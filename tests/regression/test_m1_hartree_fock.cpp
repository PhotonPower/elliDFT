// Milestone M1: restricted Hartree-Fock against the reference code.
// Acceptance criterion: total energy within 1e-8 Hartree for all test
// molecules up to benzene in STO-3G, 6-31G*, def2-SVP, cc-pVDZ, aug-cc-pVDZ.
//
// Reference: tests/reference/data/m1/hf.json (PySCF, tests/reference/gen_m1.py).
//
// The integral-direct Fock construction of M1 has no screening and no
// threads, and the development machine is not meant for long runs. Every
// test run therefore covers only the cases with at most 30 primitives (about
// one second in total). The other cases are in a test tagged [slow], run with
// scripts/build_and_test.sh --slow on a suitable machine.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "reference_data.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using nlohmann::json;

namespace {

constexpr double acceptance_threshold = 1e-8;  // Hartree
constexpr std::size_t slow_above_primitives = 30;

struct Summary {
  int cases = 0;
  double energy = 0.0, dipole = 0.0, orbital_energy = 0.0;
};

void compare_hartree_fock(bool slow_cases, Summary& summary) {
  const json reference = testing::load_reference("m1/hf.json");
  ScfOptions options;
  options.energy_tolerance = 1e-11;
  options.error_tolerance = 1e-8;

  for (const json& record : reference.at("results")) {
    const int n_functions = record.at("n_functions").get<int>();
    const std::string label = record.at("molecule").get<std::string>() + " / " +
                              record.at("basis").get<std::string>();
    const Molecule molecule = testing::molecule_from_reference(record, label);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(record.at("basis").get<std::string>()),
                                           record.at("spherical").get<bool>());
    if ((basis.n_primitives() > slow_above_primitives) != slow_cases) continue;
    INFO(label << " (" << n_functions << " functions)");
    REQUIRE(basis.n_functions() == n_functions);

    const ScfResult scf = run_rhf(molecule, basis, options);
    CHECK(scf.converged);

    // Acceptance criterion.
    const double energy_deviation = std::abs(scf.energy_total - record.at("energy_total").get<double>());
    CHECK(energy_deviation <= acceptance_threshold);

    // Further quantities; tolerances chosen here, not acceptance criteria.
    CHECK_THAT(scf.energy_nuclear_repulsion,
               WithinAbs(record.at("energy_nuclear_repulsion").get<double>(), 1e-10));
    CHECK_THAT(scf.energy_one_electron, WithinAbs(record.at("energy_one_electron").get<double>(), 1e-6));
    CHECK_THAT(scf.energy_two_electron, WithinAbs(record.at("energy_two_electron").get<double>(), 1e-6));
    const auto dipole = record.at("dipole_au").get<std::vector<double>>();
    double dipole_deviation = 0.0;
    for (std::size_t k = 0; k < 3; ++k)
      dipole_deviation = std::max(dipole_deviation, std::abs(scf.dipole[k] - dipole[k]));
    CHECK(dipole_deviation <= 1e-6);
    const auto orbital_energies = record.at("orbital_energies").get<std::vector<double>>();
    REQUIRE(orbital_energies.size() == static_cast<std::size_t>(scf.orbital_energies.size()));
    double orbital_deviation = 0.0;
    for (std::size_t i = 0; i < orbital_energies.size(); ++i)
      orbital_deviation = std::max(
          orbital_deviation, std::abs(scf.orbital_energies(static_cast<Eigen::Index>(i)) - orbital_energies[i]));
    CHECK(orbital_deviation <= 1e-6);
    CHECK(scf.n_occupied == record.at("n_electrons").get<int>() / 2);

    ++summary.cases;
    summary.energy = std::max(summary.energy, energy_deviation);
    summary.dipole = std::max(summary.dipole, dipole_deviation);
    summary.orbital_energy = std::max(summary.orbital_energy, orbital_deviation);
    // One line per case on standard output, for the milestone report.
    std::printf("HF %-8s %-12s %4d functions  %2d iterations  E = %.10f  |dE| = %.2e\n",
                record.at("molecule").get<std::string>().c_str(), record.at("basis").get<std::string>().c_str(),
                n_functions, scf.iterations, scf.energy_total, energy_deviation);
    std::fflush(stdout);
  }
}

}  // namespace

TEST_CASE("M1: Hartree-Fock energies equal the reference (up to 30 primitives)", "[m1][regression]") {
  Summary summary;
  compare_hartree_fock(false, summary);
  CHECK(summary.cases == 33);
  SUCCEED("largest deviations over " << summary.cases << " cases: energy " << summary.energy << " Hartree, dipole "
                                     << summary.dipole << ", orbital energies " << summary.orbital_energy);
}

TEST_CASE("M1: Hartree-Fock energies equal the reference (more than 30 primitives)", "[m1][regression][slow]") {
  Summary summary;
  compare_hartree_fock(true, summary);
  CHECK(summary.cases == 27);
  SUCCEED("largest deviations over " << summary.cases << " cases: energy " << summary.energy << " Hartree, dipole "
                                     << summary.dipole << ", orbital energies " << summary.orbital_energy);
}

// External anchor of the specification (M1): H2 in STO-3G at 1.4 Bohr gives
// -1.1167 Hartree (attributed there to Szabo/Ostlund; four decimals).
TEST_CASE("M1: H2 / STO-3G at 1.4 Bohr gives -1.1167 Hartree", "[m1][regression]") {
  Molecule h2;
  h2.atoms = {{1, {0.0, 0.0, 0.0}}, {1, {0.0, 0.0, 1.4}}};
  const BasisSet basis = BasisSet::build(h2, BasisLibrary::load("STO-3G"));
  const ScfResult scf = run_rhf(h2, basis);
  CHECK(scf.converged);
  CHECK(std::round(scf.energy_total * 1e4) / 1e4 == -1.1167);
  CHECK_THAT(scf.energy_nuclear_repulsion, WithinAbs(1.0 / 1.4, 1e-15));
}
