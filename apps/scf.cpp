// ellidft-scf: Hartree-Fock or Kohn-Sham calculation for a closed-shell
// molecule; prints the result as JSON.
//
// Usage: ellidft-scf <geometry.xyz> <basis name>
//                    [--functional HF|LDA|PBE|PBE0|B3LYP|B3LYP5|CAM-B3LYP]
//                    [--grid <radial>,<angular>] [--aux-basis <name>] [--guess atoms|core]
//                    [--spherical|--cartesian] [--bohr] [--charge <n>]
//
// With --aux-basis, Coulomb and exchange are built with density fitting
// (RI-J, RI-K) in that auxiliary basis, e.g. def2-universal-JKFIT.

#include <chrono>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/elements.hpp"
#include "core/parallel.hpp"
#include "core/xyz.hpp"
#include "df/df_jk.hpp"
#include "scf/rhf.hpp"

int main(int argc, char** argv) {
  using namespace ellidft;

  std::vector<std::string> positional;
  std::optional<bool> spherical;
  LengthUnit unit = LengthUnit::Angstrom;
  int charge = 0;
  std::string functional_name = "HF";
  std::string aux_name;
  std::string guess_name = "atoms";
  GridSpec grid;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--spherical")
      spherical = true;
    else if (arg == "--cartesian")
      spherical = false;
    else if (arg == "--bohr")
      unit = LengthUnit::Bohr;
    else if (arg == "--charge" && i + 1 < argc)
      charge = std::atoi(argv[++i]);
    else if (arg == "--functional" && i + 1 < argc)
      functional_name = argv[++i];
    else if (arg == "--aux-basis" && i + 1 < argc)
      aux_name = argv[++i];
    else if (arg == "--guess" && i + 1 < argc)
      guess_name = argv[++i];
    else if (arg == "--grid" && i + 1 < argc) {
      const std::string value = argv[++i];
      const auto comma = value.find(',');
      if (comma == std::string::npos) {
        std::cerr << "--grid expects <radial>,<angular>\n";
        return EXIT_FAILURE;
      }
      grid.n_radial = std::atoi(value.substr(0, comma).c_str());
      grid.n_angular = std::atoi(value.substr(comma + 1).c_str());
    }
    else if (!arg.empty() && arg[0] == '-') {
      std::cerr << "unknown option " << arg << "\n";
      return EXIT_FAILURE;
    } else
      positional.push_back(arg);
  }
  if (positional.size() != 2) {
    std::cerr << "usage: ellidft-scf <geometry.xyz> <basis name> [--functional HF|LDA|PBE|PBE0|B3LYP|B3LYP5|CAM-B3LYP] "
                 "[--grid <radial>,<angular>] [--aux-basis <name>] [--spherical|--cartesian] [--bohr] "
                 "[--charge <n>]\n";
    return EXIT_FAILURE;
  }

  try {
    Molecule molecule = read_xyz_file(positional[0], unit);
    molecule.charge = charge;
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(positional[1]), spherical);

    const auto start = std::chrono::steady_clock::now();
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
    ScfOptions options;
    if (guess_name == "core")
      options.initial_guess = InitialGuess::CoreHamiltonian;
    else if (guess_name != "atoms")
      throw std::invalid_argument("--guess expects 'atoms' or 'core'");
    const ScfResult scf = run_scf(molecule, basis, method, *jk, options);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

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
                     {"spherical", basis.spherical()},
                     {"aux_basis", aux_name},
                     {"n_functions", basis.n_functions()}};
    if (method.functional.has_density_part())
      out["method"]["grid"] = std::to_string(grid.n_radial) + "," + std::to_string(grid.n_angular) +
                              " (Treutler-Ahlrichs, Lebedev, Becke, unpruned)";
    out["scf"] = {{"converged", scf.converged},
                  {"iterations", scf.iterations},
                  {"energy_hartree", scf.energy_total},
                  {"energy_nuclear_repulsion_hartree", scf.energy_nuclear_repulsion},
                  {"energy_one_electron_hartree", scf.energy_one_electron},
                  {"energy_two_electron_hartree", scf.energy_two_electron},
                  {"energy_xc_hartree", scf.energy_xc},
                  {"electrons_on_grid", scf.electrons_on_grid},
                  {"dipole_au", {scf.dipole[0], scf.dipole[1], scf.dipole[2]}},
                  {"last_energy_change_hartree", scf.energy_change},
                  {"last_commutator_error", scf.error},
                  {"wall_time_seconds", seconds}};
    out["timings_seconds"] = {{"threads", max_threads()},
                              {"one_electron", scf.timings.one_electron},
                              {"initial_guess", scf.timings.initial_guess},
                              {"coulomb_exchange", scf.timings.coulomb_exchange},
                              {"exchange_correlation", scf.timings.exchange_correlation},
                              {"diagonalisation", scf.timings.diagonalisation},
                              {"scf_total", scf.timings.total}};
    if (const auto* df = dynamic_cast<const DensityFittedJK*>(jk.get())) {
      out["density_fitting"] = {{"auxiliary_functions", df->n_auxiliary()},
                                {"memory_gb", df->memory_bytes() / 1073741824.0},
                                {"build_seconds", df->build_seconds()}};
    }
    std::cout << out.dump(2) << "\n";
    return scf.converged ? EXIT_SUCCESS : EXIT_FAILURE;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
