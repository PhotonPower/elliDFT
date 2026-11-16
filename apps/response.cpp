// ellidft-response: ground state, polarisability on a frequency grid and
// singlet excitations of a closed-shell molecule. Prints JSON in the layout
// of the specification, section 8.
//
// Usage: ellidft-response <geometry.xyz> <basis name>
//          [--functional HF|LDA|PBE|PBE0|B3LYP|B3LYP5|CAM-B3LYP]   (default CAM-B3LYP)
//          [--grid <radial>,<angular>]                              (default 75,302)
//          [--aux-basis <name>]          density fitting for J and K
//          [--frequencies <w1,w2,...>]   photon energies in Hartree
//          [--wavelengths <l1,l2,...>]   or wavelengths in nm
//          [--gamma <g>]                 damping in Hartree (complex polarisability)
//          [--states <n>]                number of excitations (default 0)
//          [--tda]                       Tamm-Dancoff approximation instead of full response
//          [--bohr] [--charge <n>]

#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/elements.hpp"
#include "core/parallel.hpp"
#include "core/units.hpp"
#include "core/xyz.hpp"
#include "df/df_jk.hpp"
#include "response/excitations.hpp"
#include "response/polarizability.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

namespace {

std::vector<double> parse_list(const std::string& text) {
  std::vector<double> values;
  std::stringstream stream(text);
  std::string item;
  while (std::getline(stream, item, ',')) {
    std::size_t used = 0;
    const double value = std::stod(item, &used);
    if (used != item.size()) throw std::invalid_argument("invalid number '" + item + "'");
    values.push_back(value);
  }
  return values;
}

nlohmann::ordered_json matrix_json(const Eigen::Matrix3d& m) {
  nlohmann::ordered_json rows = nlohmann::ordered_json::array();
  for (int k = 0; k < 3; ++k) rows.push_back({m(k, 0), m(k, 1), m(k, 2)});
  return rows;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace ellidft;

  std::vector<std::string> positional;
  std::string functional_name = "CAM-B3LYP", aux_name;
  GridSpec grid;
  LengthUnit unit = LengthUnit::Angstrom;
  int charge = 0, n_states = 0;
  double gamma = 0.0;
  bool tda = false;
  std::vector<double> frequencies;

  try {
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      auto value = [&]() -> std::string {
        if (i + 1 >= argc) throw std::invalid_argument("option " + arg + " needs a value");
        return argv[++i];
      };
      if (arg == "--functional") functional_name = value();
      else if (arg == "--aux-basis") aux_name = value();
      else if (arg == "--grid") {
        const std::vector<double> g = parse_list(value());
        if (g.size() != 2) throw std::invalid_argument("--grid expects <radial>,<angular>");
        grid.n_radial = static_cast<int>(g[0]);
        grid.n_angular = static_cast<int>(g[1]);
      } else if (arg == "--frequencies") {
        for (double w : parse_list(value())) frequencies.push_back(w);
      } else if (arg == "--wavelengths") {
        // Photon energy in eV = 1239.841984 / wavelength in nm (specification, section 8).
        for (double nm : parse_list(value())) frequencies.push_back(1239.841984 / nm / units::hartree_in_ev);
      } else if (arg == "--gamma") gamma = std::stod(value());
      else if (arg == "--states") n_states = std::stoi(value());
      else if (arg == "--charge") charge = std::stoi(value());
      else if (arg == "--tda") tda = true;
      else if (arg == "--bohr") unit = LengthUnit::Bohr;
      else if (!arg.empty() && arg[0] == '-') throw std::invalid_argument("unknown option " + arg);
      else positional.push_back(arg);
    }
    if (positional.size() != 2) {
      std::cerr << "usage: ellidft-response <geometry.xyz> <basis name> [--functional <name>] "
                   "[--grid <radial>,<angular>] [--aux-basis <name>] [--frequencies <w,...>] "
                   "[--wavelengths <nm,...>] [--gamma <g>] [--states <n>] [--tda] [--bohr] [--charge <n>]\n";
      return EXIT_FAILURE;
    }

    Molecule molecule = read_xyz_file(positional[0], unit);
    molecule.charge = charge;
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(positional[1]));
    Method method;
    method.functional = Functional::by_name(functional_name);
    method.grid = grid;

    std::unique_ptr<JKBuilder> jk;
    if (aux_name.empty()) {
      jk = std::make_unique<DirectJK>(basis);
    } else {
      const BasisSet aux = BasisSet::build(molecule, BasisLibrary::load(aux_name));
      jk = std::make_unique<DensityFittedJK>(basis, aux);
    }

    ScfOptions scf_options;
    scf_options.energy_tolerance = 1e-10;
    scf_options.error_tolerance = 1e-8;
    const ScfResult scf = run_scf(molecule, basis, method, *jk, scf_options);
    if (!scf.converged) throw std::runtime_error("the SCF did not converge");

    const ResponseOperator op(molecule, basis, method, scf, *jk);

    nlohmann::ordered_json symbols = nlohmann::ordered_json::array();
    nlohmann::ordered_json geometry = nlohmann::ordered_json::array();
    for (const Atom& a : molecule.atoms) {
      symbols.push_back(element_symbol(a.Z));
      geometry.push_back({a.r[0], a.r[1], a.r[2]});
    }
    nlohmann::ordered_json out;
    out["molecule"] = {{"name", molecule.name}, {"charge", molecule.charge},
                       {"symbols", symbols}, {"geometry_bohr", geometry}};
    out["method"] = {{"functional", method.functional.name},
                     {"basis", basis.name()},
                     {"aux_basis", aux_name},
                     {"grid", method.functional.has_density_part()
                                  ? std::to_string(grid.n_radial) + "," + std::to_string(grid.n_angular)
                                  : std::string()},
                     {"response", tda ? "TDA" : "RPA"}};
    out["scf"] = {{"energy_hartree", scf.energy_total},
                  {"dipole_au", {scf.dipole[0], scf.dipole[1], scf.dipole[2]}}};

    const Polarizability alpha0 = polarizability(op);
    if (!alpha0.converged) throw std::runtime_error("the static response equations did not converge");
    out["polarizability_static_au"] = matrix_json(alpha0.real);

    out["polarizability_dynamic"] = nlohmann::ordered_json::array();
    for (double omega : frequencies) {
      const Polarizability alpha = polarizability(op, omega, gamma);
      if (!alpha.converged)
        throw std::runtime_error("the response equations did not converge at omega = " + std::to_string(omega));
      out["polarizability_dynamic"].push_back({{"omega_hartree", omega},
                                               {"gamma_hartree", gamma},
                                               {"alpha_re_au", matrix_json(alpha.real)},
                                               {"alpha_im_au", matrix_json(alpha.imag)}});
    }

    out["excitations"] = nlohmann::ordered_json::array();
    if (n_states > 0) {
      ExcitationOptions options;
      options.n_states = n_states;
      options.tamm_dancoff = tda;
      const ExcitationResult result = excitations(op, options);
      if (!result.converged) throw std::runtime_error("the excitation energies did not converge");
      for (const Excitation& state : result.states)
        out["excitations"].push_back(
            {{"energy_hartree", state.energy},
             {"oscillator_strength", state.oscillator_strength},
             {"transition_dipole_au",
              {state.transition_dipole[0], state.transition_dipole[1], state.transition_dipole[2]}}});
    }
    out["run"] = {{"threads", max_threads()},
                  {"scf_iterations", scf.iterations},
                  {"scf_seconds", scf.timings.total},
                  {"response_operator_applications", op.applications()}};
    std::cout << out.dump(2) << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
