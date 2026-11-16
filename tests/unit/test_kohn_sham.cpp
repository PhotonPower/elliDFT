// Kohn-Sham: checks that need no reference code (specification, section 7).
// Small basis sets and grids keep these tests short.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "basis/basis_set.hpp"
#include "core/parallel.hpp"
#include "core/xyz.hpp"
#include "scf/rhf.hpp"
#include "scf/xc_matrix.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Molecule load(const std::string& stem) {
  return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/" + stem + ".xyz");
}

// A physically reasonable density matrix: Hartree-Fock in the same basis.
Eigen::MatrixXd some_density(const Molecule& molecule, const BasisSet& basis) {
  return run_rhf(molecule, basis).density;
}

Eigen::MatrixXd symmetric_perturbation(Eigen::Index n) {
  Eigen::MatrixXd m(n, n);
  unsigned state = 12345u;
  for (Eigen::Index i = 0; i < n; ++i)
    for (Eigen::Index j = 0; j <= i; ++j) {
      state = state * 1664525u + 1013904223u;
      m(i, j) = m(j, i) = static_cast<double>(state >> 8) / static_cast<double>(1u << 24) - 0.5;
    }
  return m;
}

Eigen::Matrix3d rotation(const Eigen::Vector3d& axis, double angle) {
  const Eigen::Vector3d u = axis.normalized();
  Eigen::Matrix3d cross;
  cross << 0.0, -u.z(), u.y(), u.z(), 0.0, -u.x(), -u.y(), u.x(), 0.0;
  return Eigen::Matrix3d::Identity() + std::sin(angle) * cross + (1.0 - std::cos(angle)) * cross * cross;
}

}  // namespace

// Exchange-correlation matrix against the numerical derivative of the
// exchange-correlation energy with respect to the density matrix:
// dE_xc = sum V_mu,nu dD_mu,nu.
TEST_CASE("XC matrix is the derivative of the XC energy", "[scf][xc_matrix]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("6-31G*"));
  const MolecularGrid grid = build_molecular_grid(water, {30, 110});
  const Eigen::MatrixXd D = some_density(water, basis);
  const Eigen::MatrixXd dD = symmetric_perturbation(basis.n_functions());

  for (const char* name : {"LDA", "PBE"}) {
    INFO(name);
    const auto xc = make_libxc_evaluator(Functional::by_name(name));
    const XcIntegrals at_D = integrate_xc(basis, grid, *xc, D);
    CHECK((at_D.matrix - at_D.matrix.transpose()).cwiseAbs().maxCoeff() < 1e-14);
    const double analytic = at_D.matrix.cwiseProduct(dD).sum();

    auto difference = [&](double h) {
      return (integrate_xc(basis, grid, *xc, D + h * dD).energy - integrate_xc(basis, grid, *xc, D - h * dD).energy) /
             (2.0 * h);
    };
    const double h = 1e-4;
    const double numeric = (4.0 * difference(h) - difference(2.0 * h)) / 3.0;
    CHECK_THAT(numeric, WithinRel(analytic, 1e-7));
  }
}

TEST_CASE("XC integration: own LDA, block size, unsymmetric density", "[scf][xc_matrix]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("STO-3G"));
  const MolecularGrid grid = build_molecular_grid(water, {30, 110});
  const Eigen::MatrixXd D = some_density(water, basis);
  const auto libxc_lda = make_libxc_evaluator(Functional::by_name("LDA"));
  const auto own_lda = make_own_lda_evaluator();

  const XcIntegrals a = integrate_xc(basis, grid, *libxc_lda, D);
  const XcIntegrals b = integrate_xc(basis, grid, *own_lda, D);
  CHECK_THAT(b.energy, WithinAbs(a.energy, 1e-11));
  CHECK((b.matrix - a.matrix).cwiseAbs().maxCoeff() < 1e-11);

  // The result does not depend on the block size beyond rounding.
  const XcIntegrals c = integrate_xc(basis, grid, *libxc_lda, D, 37);
  CHECK_THAT(c.energy, WithinAbs(a.energy, 1e-12));
  CHECK((c.matrix - a.matrix).cwiseAbs().maxCoeff() < 1e-12);

  // Only the symmetric part of the density enters.
  Eigen::MatrixXd skew = symmetric_perturbation(basis.n_functions());
  skew = skew.triangularView<Eigen::StrictlyUpper>().toDenseMatrix() -
         skew.triangularView<Eigen::StrictlyUpper>().toDenseMatrix().transpose();
  const auto pbe = make_libxc_evaluator(Functional::by_name("PBE"));
  const XcIntegrals d = integrate_xc(basis, grid, *pbe, D);
  const XcIntegrals e = integrate_xc(basis, grid, *pbe, D + skew);
  CHECK_THAT(e.energy, WithinAbs(d.energy, 1e-12));
  CHECK((e.matrix - d.matrix).cwiseAbs().maxCoeff() < 1e-12);

  CHECK_THROWS_AS(integrate_xc(basis, grid, *pbe, Eigen::MatrixXd::Zero(2, 2)), std::invalid_argument);
  CHECK_THROWS_AS(integrate_xc(basis, grid, *pbe, D, 0), std::invalid_argument);
}

// Acceptance criterion M2: integrated number of electrons within 1e-6 per
// electron on the standard grid (75, 302); the (99, 590) grid is checked in
// the regression test. (Relative criterion by decision of the client of
// 2026-11-14, see docs/Entscheidungen.md.)
TEST_CASE("Kohn-Sham: electron count on the standard grid and density conditions", "[scf][ks][m2]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("6-31G*"));
  for (const char* name : {"LDA", "PBE"}) {
    INFO(name);
    Method method;
    method.functional = Functional::by_name(name);
    method.grid = {75, 302};
    const ScfResult scf = run_scf(water, basis, method);
    REQUIRE(scf.converged);
    CHECK(std::abs(scf.electrons_on_grid - 10.0) / 10.0 <= 1e-6);
    const Eigen::MatrixXd& D = scf.density;
    const Eigen::MatrixXd& S = scf.overlap;
    CHECK_THAT((D * S).trace(), WithinAbs(10.0, 1e-10));
    CHECK((D * S * D - 2.0 * D).cwiseAbs().maxCoeff() < 1e-10);
    CHECK((scf.fock * D * S - S * D * scf.fock).cwiseAbs().maxCoeff() < 1e-7);
    CHECK(scf.energy_xc < 0.0);
    CHECK_THAT(scf.energy_total,
               WithinAbs(scf.energy_one_electron + scf.energy_two_electron + scf.energy_xc +
                             scf.energy_nuclear_repulsion, 1e-11));
  }
}

// Acceptance criterion M2: an arbitrary rigid rotation changes the energy by
// at most 1e-6 Hartree. (The grid is not rotationally invariant.)
TEST_CASE("Kohn-Sham: rotational invariance of the energy", "[scf][ks][m2]") {
  const Molecule water = load("h2o");
  const Eigen::Matrix3d R = rotation({0.4, -0.9, 0.6}, 0.77);
  Molecule rotated = water;
  for (Atom& atom : rotated.atoms) {
    const Eigen::Vector3d r = R * Eigen::Vector3d(atom.r[0], atom.r[1], atom.r[2]) + Eigen::Vector3d(0.3, 1.1, -0.6);
    atom.r = {r.x(), r.y(), r.z()};
  }
  const BasisLibrary library = BasisLibrary::load("6-31G*");
  for (const char* name : {"LDA", "PBE"}) {
    INFO(name);
    Method method;
    method.functional = Functional::by_name(name);
    method.grid = {75, 302};
    const ScfResult a = run_scf(water, BasisSet::build(water, library), method);
    const ScfResult b = run_scf(rotated, BasisSet::build(rotated, library), method);
    REQUIRE(a.converged);
    REQUIRE(b.converged);
    INFO("energy change " << b.energy_total - a.energy_total);
    CHECK(std::abs(b.energy_total - a.energy_total) <= 1e-6);
    SUCCEED(name << ": energy change under rotation " << b.energy_total - a.energy_total);
  }
}

// Rule 11 of the specification: the energy depends on the number of threads
// by less than 1e-10 Hartree.
TEST_CASE("energies do not depend on the number of threads", "[scf][ks][threads]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("6-31G*"));
  const int available = max_threads();
  for (const char* name : {"HF", "PBE"}) {
    INFO(name);
    Method method;
    method.functional = Functional::by_name(name);
    method.grid = {40, 110};
    set_num_threads(1);
    const ScfResult serial = run_scf(water, basis, method);
    set_num_threads(std::max(available, 4));
    const ScfResult threaded = run_scf(water, basis, method);
    set_num_threads(available);
    REQUIRE(serial.converged);
    REQUIRE(threaded.converged);
    INFO("difference " << threaded.energy_total - serial.energy_total);
    CHECK(std::abs(threaded.energy_total - serial.energy_total) < 1e-10);
  }
}

TEST_CASE("run_scf: Hartree-Fock limit and invalid functionals", "[scf][ks]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("STO-3G"));
  const ScfResult hf = run_scf(water, basis, Method{});
  CHECK_THAT(hf.energy_total, WithinAbs(run_rhf(water, basis).energy_total, 1e-12));
  CHECK(hf.energy_xc == 0.0);

  // Long-range exchange without a range separation parameter is rejected.
  Method range_separated;
  range_separated.functional = Functional::by_name("PBE");
  range_separated.functional.long_range_exchange = 0.46;
  range_separated.functional.omega = 0.0;
  CHECK_THROWS_AS(run_scf(water, basis, range_separated), std::invalid_argument);
}
