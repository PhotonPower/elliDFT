// Milestone M2: Kohn-Sham total energies against the reference code.
// Acceptance criteria:
//   - LDA and PBE total energies within 5e-6 Hartree, both programs on the
//     unpruned (99, 590) grid;
//   - integrated number of electrons within 1e-7 per electron on the
//     (99, 590) grid. (Changed by the client on 2026-11-14 from an absolute
//     error of 1e-6, which grows with the size of the molecule and is not met
//     by the reference code for benzene; see docs/Entscheidungen.md.)
//
// Reference: tests/reference/data/m2/ks.json (PySCF, tests/reference/gen_m2.py).
//
// Every test run covers the atoms and H2 (about one second). The other cases
// are tagged [slow] and meant for a server (scripts/build_and_test.sh --slow).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "reference_data.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using nlohmann::json;

namespace {

constexpr double energy_threshold = 5e-6;     // Hartree
constexpr double electron_threshold_per_electron = 1e-7;

struct Summary {
  int cases = 0;
  double energy = 0.0, electrons = 0.0, xc = 0.0, dipole = 0.0;
};

bool is_quick(const BasisSet& basis, const Molecule& molecule) {
  return molecule.n_atoms() <= 2 && basis.n_primitives() <= 12;
}

void compare_kohn_sham(bool slow_cases, Summary& summary) {
  const json reference = testing::load_reference("m2/ks.json");
  const json& grid = reference.at("meta").at("settings").at("grid");
  Method method;
  method.grid = {grid.at("n_radial").get<int>(), grid.at("n_angular").get<int>()};
  REQUIRE(method.grid.n_radial == 99);
  REQUIRE(method.grid.n_angular == 590);
  ScfOptions options;
  options.energy_tolerance = 1e-10;
  options.error_tolerance = 1e-7;

  for (const json& record : reference.at("results")) {
    const std::string label = record.at("molecule").get<std::string>() + " / " +
                              record.at("basis").get<std::string>() + " / " +
                              record.at("functional").get<std::string>();
    const Molecule molecule = testing::molecule_from_reference(record, label);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(record.at("basis").get<std::string>()),
                                           record.at("spherical").get<bool>());
    if (is_quick(basis, molecule) == slow_cases) continue;
    INFO(label);
    REQUIRE(basis.n_functions() == record.at("n_functions").get<int>());

    method.functional = Functional::by_name(record.at("functional").get<std::string>());
    const ScfResult scf = run_scf(molecule, basis, method, options);
    CHECK(scf.converged);

    const double energy_deviation = std::abs(scf.energy_total - record.at("energy_total").get<double>());
    const int n_electrons = record.at("n_electrons").get<int>();
    const double electron_error = std::abs(scf.electrons_on_grid - n_electrons) / n_electrons;
    CHECK(energy_deviation <= energy_threshold);                 // acceptance criterion
    CHECK(electron_error <= electron_threshold_per_electron);    // acceptance criterion

    // Further comparisons; tolerances chosen here.
    const double xc_deviation = std::abs(scf.energy_xc - record.at("energy_xc").get<double>());
    CHECK(xc_deviation <= 1e-6);
    CHECK(std::abs(scf.electrons_on_grid - record.at("electrons_on_grid").get<double>()) <= 1e-8);
    const auto dipole = record.at("dipole_au").get<std::vector<double>>();
    double dipole_deviation = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
      dipole_deviation = std::max(dipole_deviation, std::abs(scf.dipole[k] - dipole[k]));
    CHECK(dipole_deviation <= 1e-5);

    ++summary.cases;
    summary.energy = std::max(summary.energy, energy_deviation);
    summary.electrons = std::max(summary.electrons, electron_error);
    summary.xc = std::max(summary.xc, xc_deviation);
    summary.dipole = std::max(summary.dipole, dipole_deviation);
    std::printf("KS %-8s %-12s %-4s %2d iterations  E = %.10f  |dE| = %.2e  |N(grid) - N| / N = %.2e\n",
                record.at("molecule").get<std::string>().c_str(), record.at("basis").get<std::string>().c_str(),
                record.at("functional").get<std::string>().c_str(), scf.iterations, scf.energy_total,
                energy_deviation, electron_error);
    std::fflush(stdout);
  }
}

}  // namespace

TEST_CASE("M2: LDA and PBE energies equal the reference (atoms and H2)", "[m2][regression]") {
  Summary summary;
  compare_kohn_sham(false, summary);
  CHECK(summary.cases == 6);
  SUCCEED("largest deviations over " << summary.cases << " cases: energy " << summary.energy
                                     << " Hartree, electrons (relative) " << summary.electrons << ", E_xc " << summary.xc
                                     << ", dipole " << summary.dipole);
}

TEST_CASE("M2: LDA and PBE energies equal the reference (all other cases)", "[m2][regression][slow]") {
  Summary summary;
  compare_kohn_sham(true, summary);
  CHECK(summary.cases == 26);
  SUCCEED("largest deviations over " << summary.cases << " cases: energy " << summary.energy
                                     << " Hartree, electrons (relative) " << summary.electrons << ", E_xc " << summary.xc
                                     << ", dipole " << summary.dipole);
}
