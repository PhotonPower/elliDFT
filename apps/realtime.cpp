// ellidft-realtime: polarisability of a closed-shell molecule from the
// propagation of its orbitals in time after a short field pulse.
// Prints JSON; "polarizability_dynamic" has the layout of ellidft-response.
//
// Usage: ellidft-realtime <geometry.xyz> <basis name>
//          [--functional HF|LDA|PBE|PBE0|B3LYP|B3LYP5|CAM-B3LYP]   (default CAM-B3LYP)
//          [--grid <radial>,<angular>]                              (default 75,302)
//          [--aux-basis <name>]          density fitting for J and K
//          [--time-step <dt>]            atomic units of time (default 0.1)
//          [--time <T>]                  length of the propagation; default 12 / gamma
//          [--kick <kappa>]              strength of the pulse (default 1e-4)
//          [--one-sided]                 one propagation per direction instead of +kappa and -kappa
//          [--directions <xyz>]          pulse directions (default xyz)
//          [--propagator magnus4|composition4|midpoint]
//          [--gamma <g>]                 damping in Hartree (required, > 0)
//          [--frequencies <w1,w2,...>]   photon energies in Hartree
//          [--wavelengths <l1,l2,...>]   or wavelengths in nm
//          [--spectrum <from>,<to>,<n>]  or n equidistant photon energies in Hartree
//          [--compare]                   also solve the linear response equations and report the difference
//          [--signal <file>]             write d mu_i / d kappa_j (t) as text
//          [--bohr] [--charge <n>]

#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
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
#include "realtime/propagator.hpp"
#include "realtime/spectrum.hpp"
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

template <typename Matrix>
nlohmann::ordered_json matrix_json(const Matrix& m) {
  nlohmann::ordered_json rows = nlohmann::ordered_json::array();
  for (int k = 0; k < 3; ++k) rows.push_back({m(k, 0), m(k, 1), m(k, 2)});
  return rows;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace ellidft;

  std::vector<std::string> positional;
  std::string functional_name = "CAM-B3LYP", aux_name, signal_file, directions = "xyz", propagator_name = "magnus4";
  GridSpec grid;
  LengthUnit unit = LengthUnit::Angstrom;
  int charge = 0;
  double gamma = 0.0, total_time = 0.0;
  bool compare = false;
  RealTimeOptions options;
  KickOptions kick;
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
      } else if (arg == "--spectrum") {
        const std::vector<double> s = parse_list(value());
        if (s.size() != 3 || s[2] < 2.0 || !(s[1] > s[0])) throw std::invalid_argument("--spectrum expects <from>,<to>,<n>");
        const int n = static_cast<int>(s[2]);
        for (int k = 0; k < n; ++k) frequencies.push_back(s[0] + (s[1] - s[0]) * k / (n - 1));
      } else if (arg == "--gamma") gamma = std::stod(value());
      else if (arg == "--time-step") options.time_step = std::stod(value());
      else if (arg == "--time") total_time = std::stod(value());
      else if (arg == "--kick") kick.strength = std::stod(value());
      else if (arg == "--one-sided") kick.two_sided = false;
      else if (arg == "--directions") directions = value();
      else if (arg == "--propagator") propagator_name = value();
      else if (arg == "--compare") compare = true;
      else if (arg == "--signal") signal_file = value();
      else if (arg == "--charge") charge = std::stoi(value());
      else if (arg == "--bohr") unit = LengthUnit::Bohr;
      else if (!arg.empty() && arg[0] == '-') throw std::invalid_argument("unknown option " + arg);
      else positional.push_back(arg);
    }
    if (positional.size() != 2) {
      std::cerr << "usage: ellidft-realtime <geometry.xyz> <basis name> --gamma <g> [--functional <name>] "
                   "[--grid <radial>,<angular>] [--aux-basis <name>] [--time-step <dt>] [--time <T>] [--kick <kappa>] "
                   "[--one-sided] [--directions <xyz>] [--propagator magnus4|composition4|midpoint] [--frequencies <w,...>] "
                   "[--wavelengths <nm,...>] [--spectrum <from>,<to>,<n>] [--compare] [--signal <file>] [--bohr] "
                   "[--charge <n>]\n";
      return EXIT_FAILURE;
    }
    if (!(gamma > 0.0)) throw std::invalid_argument("--gamma is required and must be positive");
    if (!(options.time_step > 0.0)) throw std::invalid_argument("--time-step must be positive");
    if (propagator_name == "magnus4") options.propagator = Propagator::Magnus4;
    else if (propagator_name == "midpoint") options.propagator = Propagator::Midpoint;
    else if (propagator_name == "composition4") options.propagator = Propagator::Composition4;
    else throw std::invalid_argument("--propagator expects 'magnus4', 'composition4' or 'midpoint'");
    kick.directions = {false, false, false};
    for (char c : directions) {
      if (c != 'x' && c != 'y' && c != 'z') throw std::invalid_argument("--directions expects letters out of xyz");
      kick.directions[static_cast<std::size_t>(c - 'x')] = true;
    }
    if (total_time == 0.0) total_time = 12.0 / gamma;   // neglected tail exp(-12) = 6e-6
    if (!(total_time > 0.0)) throw std::invalid_argument("--time must be positive");
    options.n_steps = static_cast<int>(std::ceil(total_time / options.time_step));

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

    // The propagation starts from the ground state; what is left of the SCF
    // error moves in time like a response. Hence tighter than elsewhere.
    ScfOptions scf_options;
    scf_options.energy_tolerance = 1e-11;
    scf_options.error_tolerance = 1e-9;
    const ScfResult scf = run_scf(molecule, basis, method, *jk, scf_options);
    if (!scf.converged) throw std::runtime_error("the SCF did not converge");

    const RealTimePropagator propagator(molecule, basis, method, scf, *jk);
    std::cerr << "real-time propagation: " << options.n_steps << " steps of " << options.time_step << " a.u., "
              << (kick.two_sided ? 2 : 1) * static_cast<int>(directions.size()) << " propagation(s)\n";
    const DipoleResponse response = dipole_response(propagator, options, kick);
    const double tail = std::exp(-gamma * response.total_time());
    std::cerr << "  " << response.fock_builds << " Fock matrices, energy drift " << std::scientific
              << std::setprecision(1) << response.energy_drift << " Hartree, neglected tail of the signal " << tail
              << "\n";
    if (!response.steps_converged) std::cerr << "  WARNING: the iteration of some time steps did not converge\n";
    if (tail > 1e-3) std::cerr << "  WARNING: the propagation is short for this damping (tail above 1e-3)\n";

    if (!signal_file.empty()) {
      std::ofstream signal(signal_file);
      if (!signal) throw std::runtime_error("cannot write " + signal_file);
      signal << "# d mu_i / d kappa_j in atomic units; pulse directions " << directions << "\n"
             << "# t_au xx xy xz yx yy yz zx zy zz\n"
             << std::scientific << std::setprecision(10);
      for (std::size_t n = 0; n < response.response.size(); ++n) {
        signal << response.time_step * static_cast<double>(n);
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) signal << " " << response.response[n](a, b);
        signal << "\n";
      }
    }

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
                     {"response", "real-time propagation"}};
    out["scf"] = {{"energy_hartree", scf.energy_total},
                  {"dipole_au", {scf.dipole[0], scf.dipole[1], scf.dipole[2]}}};
    out["realtime"] = {{"propagator", propagator_name},
                       {"time_step_au", options.time_step},
                       {"steps", options.n_steps},
                       {"total_time_au", response.total_time()},
                       {"kick_au", kick.strength},
                       {"two_sided", kick.two_sided},
                       {"directions", directions},
                       {"propagations", response.propagations},
                       {"fock_builds", response.fock_builds},
                       {"steps_converged", response.steps_converged},
                       {"energy_drift_hartree", response.energy_drift},
                       {"orthonormality_error", response.orthonormality_error},
                       {"neglected_tail", tail}};

    std::unique_ptr<ResponseOperator> op;
    if (compare) op = std::make_unique<ResponseOperator>(molecule, basis, method, scf, *jk);
    double largest_deviation = 0.0;
    out["polarizability_dynamic"] = nlohmann::ordered_json::array();
    for (double omega : frequencies) {
      const Eigen::Matrix3cd alpha = polarizability(response, omega, gamma);
      nlohmann::ordered_json entry = {{"omega_hartree", omega},
                                      {"gamma_hartree", gamma},
                                      {"alpha_re_au", matrix_json(alpha.real())},
                                      {"alpha_im_au", matrix_json(alpha.imag())}};
      if (compare) {
        const Polarizability reference = polarizability(*op, omega, gamma);
        if (!reference.converged)
          throw std::runtime_error("the response equations did not converge at omega = " + std::to_string(omega));
        double scale = 0.0, difference = 0.0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) {
            if (!kick.directions[static_cast<std::size_t>(b)]) continue;
            const std::complex<double> ref(reference.real(a, b), reference.imag(a, b));
            scale = std::max(scale, std::abs(ref));
            difference = std::max(difference, std::abs(alpha(a, b) - ref));
          }
        entry["linear_response_re_au"] = matrix_json(reference.real);
        entry["linear_response_im_au"] = matrix_json(reference.imag);
        entry["relative_deviation"] = difference / scale;
        largest_deviation = std::max(largest_deviation, difference / scale);
      }
      out["polarizability_dynamic"].push_back(entry);
    }
    if (compare && !frequencies.empty()) {
      out["realtime"]["largest_relative_deviation_from_linear_response"] = largest_deviation;
      std::cerr << "  agreement with linear response: " << largest_deviation << " (largest relative deviation)\n";
    }
    out["run"] = {{"threads", max_threads()}, {"scf_iterations", scf.iterations}, {"scf_seconds", scf.timings.total}};
    std::cout << out.dump(2) << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
