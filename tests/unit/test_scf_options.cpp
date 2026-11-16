// SCF: start density, linear dependencies, density differences, timings.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "core/xyz.hpp"
#include "integrals/one_electron.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;

namespace {

Molecule load(const std::string& stem) {
  return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/" + stem + ".xyz");
}

}  // namespace

TEST_CASE("superposition of atomic densities", "[scf][guess]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("6-31G*"));
  const Eigen::MatrixXd D = atomic_density_guess(water, basis);
  const Eigen::MatrixXd S = overlap_matrix(basis);
  const Eigen::Index n = basis.n_functions();
  REQUIRE(D.rows() == n);

  CHECK((D - D.transpose()).cwiseAbs().maxCoeff() < 1e-12);
  // Blocks between different atoms are zero; the atomic blocks hold Z electrons.
  // 6-31G* (cartesian): O has 15 functions, each H has 2.
  CHECK(D.block(0, 15, 15, 4).cwiseAbs().maxCoeff() == 0.0);
  CHECK(D.block(15, 17, 2, 2).cwiseAbs().maxCoeff() == 0.0);
  CHECK_THAT((D.block(0, 0, 15, 15) * S.block(0, 0, 15, 15)).trace(), WithinAbs(8.0, 1e-8));
  CHECK_THAT((D.block(15, 15, 2, 2) * S.block(15, 15, 2, 2)).trace(), WithinAbs(1.0, 1e-8));
  // Both hydrogen atoms get the same density.
  CHECK((D.block(15, 15, 2, 2) - D.block(17, 17, 2, 2)).cwiseAbs().maxCoeff() < 1e-12);
  // Spherical average: the three p directions of oxygen are equivalent.
  // Oxygen in 6-31G*: three s functions (0..2), then px, py, pz of the first
  // p shell (3..5).
  CHECK_THAT(D(3, 3), WithinAbs(D(4, 4), 1e-8));
  CHECK_THAT(D(4, 4), WithinAbs(D(5, 5), 1e-8));
  CHECK_THAT(D(3, 4), WithinAbs(0.0, 1e-8));
  CHECK_THAT(D(0, 3), WithinAbs(0.0, 1e-8));  // no s-p mixing in a spherical atom
}

TEST_CASE("closed-shell atom: the atomic density is the converged density", "[scf][guess]") {
  Molecule neon;
  neon.atoms = {{10, {0.3, -0.1, 0.2}}};
  const BasisSet basis = BasisSet::build(neon, BasisLibrary::load("def2-SVP"));
  ScfOptions core;
  core.initial_guess = InitialGuess::CoreHamiltonian;
  const ScfResult from_core = run_rhf(neon, basis, core);
  const ScfResult from_atoms = run_rhf(neon, basis);
  REQUIRE(from_core.converged);
  REQUIRE(from_atoms.converged);
  CHECK_THAT(from_atoms.energy_total, WithinAbs(from_core.energy_total, 1e-9));
  CHECK(from_atoms.iterations <= 3);
  CHECK(from_atoms.iterations < from_core.iterations);
}

TEST_CASE("start density does not change the result", "[scf][guess]") {
  for (const char* stem : {"h2o", "hcl"}) {
    INFO(stem);
    const Molecule molecule = load(stem);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("6-31G*"));
    ScfOptions core;
    core.initial_guess = InitialGuess::CoreHamiltonian;
    for (const char* name : {"HF", "PBE"}) {
      INFO(name);
      Method method;
      method.functional = Functional::by_name(name);
      method.grid = {40, 110};
      const ScfResult a = run_scf(molecule, basis, method, core);
      const ScfResult b = run_scf(molecule, basis, method);
      REQUIRE(a.converged);
      REQUIRE(b.converged);
      INFO("iterations: core " << a.iterations << ", atomic densities " << b.iterations);
      CHECK_THAT(b.energy_total, WithinAbs(a.energy_total, 1e-8));
      CHECK(b.iterations <= a.iterations);
    }
  }
}

TEST_CASE("density differences do not change the result", "[scf]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("6-31G*"));
  ScfOptions full, incremental, never_rebuilt;
  full.difference_densities = false;
  never_rebuilt.full_rebuild_interval = 0;
  const ScfResult a = run_rhf(water, basis, full);
  const ScfResult b = run_rhf(water, basis, incremental);
  const ScfResult c = run_rhf(water, basis, never_rebuilt);
  REQUIRE(a.converged);
  REQUIRE(b.converged);
  REQUIRE(c.converged);
  CHECK_THAT(b.energy_total, WithinAbs(a.energy_total, 1e-9));
  CHECK_THAT(c.energy_total, WithinAbs(a.energy_total, 1e-9));
  CHECK(b.iterations == a.iterations);
}

TEST_CASE("near-linear dependencies are removed from the orbital space", "[scf][lindep]") {
  // Helium with two s functions, and the same with a third one that almost
  // duplicates the second.
  auto shell = [](double exponent) {
    Shell s;
    s.atom = 0;
    s.center = {0.0, 0.0, 0.0};
    s.l = 0;
    s.exponents = {exponent};
    s.coefficients = {1.0};
    return s;
  };
  Molecule helium;
  helium.atoms = {{2, {0.0, 0.0, 0.0}}};
  const BasisSet two = BasisSet::from_shells("two", true, {shell(3.0), shell(0.6)});
  const BasisSet three = BasisSet::from_shells("three", true, {shell(3.0), shell(0.6), shell(0.6 * (1.0 + 1e-6))});

  const ScfResult reference = run_rhf(helium, two);
  REQUIRE(reference.converged);
  CHECK(reference.n_orbitals == 2);

  const ScfResult reduced = run_rhf(helium, three);
  REQUIRE(reduced.converged);
  CHECK(reduced.overlap_min_eigenvalue < 1e-8);
  CHECK(reduced.n_orbitals == 2);          // one combination removed
  CHECK(reduced.orbital_energies.size() == 2);
  CHECK(reduced.coefficients.rows() == 3);
  CHECK(reduced.coefficients.cols() == 2);
  CHECK_THAT(reduced.energy_total, WithinAbs(reference.energy_total, 1e-6));

  // A threshold of zero keeps every function.
  ScfOptions keep_all;
  keep_all.overlap_threshold = 0.0;
  keep_all.max_iterations = 3;
  CHECK(run_rhf(helium, three, keep_all).n_orbitals == 3);

  // More occupied orbitals than independent functions.
  Molecule beryllium;
  beryllium.atoms = {{4, {0.0, 0.0, 0.0}}};
  const BasisSet duplicate = BasisSet::from_shells("dup", true, {shell(1.0), shell(1.0 * (1.0 + 1e-9))});
  CHECK_THROWS_AS(run_rhf(beryllium, duplicate), std::invalid_argument);
}

TEST_CASE("SCF reports the time of its parts", "[scf][timings]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("STO-3G"));
  Method method;
  method.functional = Functional::by_name("PBE");
  method.grid = {20, 50};
  const ScfResult scf = run_scf(water, basis, method);
  const ScfTimings& t = scf.timings;
  CHECK(t.one_electron > 0.0);
  CHECK(t.initial_guess > 0.0);
  CHECK(t.coulomb_exchange > 0.0);
  CHECK(t.exchange_correlation > 0.0);
  CHECK(t.diagonalisation > 0.0);
  CHECK(t.total >= t.one_electron + t.initial_guess + t.coulomb_exchange + t.exchange_correlation);
  // Hartree-Fock has no exchange-correlation part.
  CHECK(run_rhf(water, basis).timings.exchange_correlation == 0.0);
}
