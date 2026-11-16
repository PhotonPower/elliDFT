// Restricted Hartree-Fock: checks that need no reference code
// (specification, section 7, "Proben ohne Referenz").

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "basis/basis_set.hpp"
#include "core/xyz.hpp"
#include "scf/diis.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;

namespace {

Molecule load(const std::string& stem) {
  return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/" + stem + ".xyz");
}

ScfOptions tight() {
  ScfOptions options;
  options.energy_tolerance = 1e-12;
  options.error_tolerance = 1e-9;
  return options;
}

// Proper rotation about an arbitrary axis (Rodrigues' formula, written out).
Eigen::Matrix3d rotation(const Eigen::Vector3d& axis, double angle) {
  const Eigen::Vector3d u = axis.normalized();
  Eigen::Matrix3d cross;
  cross << 0.0, -u.z(), u.y(), u.z(), 0.0, -u.x(), -u.y(), u.x(), 0.0;
  return Eigen::Matrix3d::Identity() + std::sin(angle) * cross + (1.0 - std::cos(angle)) * cross * cross;
}

Molecule transformed(const Molecule& molecule, const Eigen::Matrix3d& R, const Eigen::Vector3d& shift) {
  Molecule out = molecule;
  for (Atom& atom : out.atoms) {
    const Eigen::Vector3d r = R * Eigen::Vector3d(atom.r[0], atom.r[1], atom.r[2]) + shift;
    atom.r = {r.x(), r.y(), r.z()};
  }
  return out;
}

}  // namespace

TEST_CASE("RHF: density matrix and stationarity conditions", "[scf][rhf]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("def2-SVP"));
  const ScfResult scf = run_rhf(water, basis, tight());
  REQUIRE(scf.converged);
  const Eigen::MatrixXd& D = scf.density;
  const Eigen::MatrixXd& S = scf.overlap;
  const Eigen::MatrixXd& F = scf.fock;
  const Eigen::MatrixXd& C = scf.coefficients;

  // Tr(DS) = number of electrons, D S D = 2 D.
  CHECK_THAT((D * S).trace(), WithinAbs(10.0, 1e-11));
  CHECK((D * S * D - 2.0 * D).cwiseAbs().maxCoeff() < 1e-11);
  CHECK((D - D.transpose()).cwiseAbs().maxCoeff() < 1e-13);
  // Orbitals are orthonormal and diagonalise the Fock matrix.
  const Eigen::Index n = basis.n_functions();
  CHECK((C.transpose() * S * C - Eigen::MatrixXd::Identity(n, n)).cwiseAbs().maxCoeff() < 1e-11);
  const Eigen::MatrixXd F_mo = C.transpose() * F * C;
  CHECK((F_mo - Eigen::MatrixXd(scf.orbital_energies.asDiagonal())).cwiseAbs().maxCoeff() < 1e-10);
  // Stationarity: F D S = S D F.
  CHECK((F * D * S - S * D * F).cwiseAbs().maxCoeff() < 1e-8);
  // Aufbau: orbital energies ascending, gap between occupied and virtual.
  REQUIRE(scf.n_occupied == 5);
  for (Eigen::Index i = 1; i < n; ++i) CHECK(scf.orbital_energies(i) >= scf.orbital_energies(i - 1));
  CHECK(scf.orbital_energies(4) < 0.0);
  CHECK(scf.orbital_energies(5) > scf.orbital_energies(4));
  // Energy bookkeeping.
  CHECK_THAT(scf.energy_total,
             WithinAbs(scf.energy_one_electron + scf.energy_two_electron + scf.energy_nuclear_repulsion, 1e-12));
  CHECK_THAT(scf.energy_nuclear_repulsion, WithinAbs(water.nuclear_repulsion_energy(), 0.0));
  CHECK(scf.overlap_min_eigenvalue > 0.0);
  // Water in this orientation: dipole along z only (C2 axis).
  CHECK_THAT(scf.dipole[0], WithinAbs(0.0, 1e-9));
  CHECK_THAT(scf.dipole[1], WithinAbs(0.0, 1e-9));
  CHECK(std::abs(scf.dipole[2]) > 0.5);
}

TEST_CASE("RHF: invariance under rotation and translation", "[scf][rhf]") {
  const Molecule water_molecule = load("h2o");  // small, to keep the test short
  const Eigen::Matrix3d R = rotation({0.3, -1.1, 0.7}, 1.234);
  const Eigen::Vector3d shift(1.9, -0.4, 2.7);
  REQUIRE((R * R.transpose() - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() < 1e-14);
  const Molecule moved = transformed(water_molecule, R, shift);

  for (const char* basis_name : {"def2-SVP", "6-31G*"}) {  // spherical and cartesian d functions
    INFO(basis_name);
    const BasisLibrary library = BasisLibrary::load(basis_name);
    const ScfResult a = run_rhf(water_molecule, BasisSet::build(water_molecule, library), tight());
    const ScfResult b = run_rhf(moved, BasisSet::build(moved, library), tight());
    REQUIRE(a.converged);
    REQUIRE(b.converged);
    CHECK_THAT(b.energy_total, WithinAbs(a.energy_total, 1e-9));
    CHECK_THAT(b.energy_one_electron, WithinAbs(a.energy_one_electron, 1e-7));
    // The dipole of a neutral molecule rotates with it and ignores the shift.
    const Eigen::Vector3d rotated = R * Eigen::Vector3d(a.dipole[0], a.dipole[1], a.dipole[2]);
    for (int k = 0; k < 3; ++k) CHECK_THAT(b.dipole[static_cast<std::size_t>(k)], WithinAbs(rotated(k), 1e-7));
    for (Eigen::Index i = 0; i < a.orbital_energies.size(); ++i)
      CHECK_THAT(b.orbital_energies(i), WithinAbs(a.orbital_energies(i), 1e-7));
  }
}

TEST_CASE("RHF: size consistency of two distant molecules", "[scf][rhf]") {
  const Molecule hf = load("hf"), h2 = load("h2");
  // H2 moved 200 Bohr away and turned.
  const Molecule far_h2 = transformed(h2, rotation({1.0, 2.0, 0.5}, 0.8), {200.0, 30.0, -10.0});
  Molecule both = hf;
  both.atoms.insert(both.atoms.end(), far_h2.atoms.begin(), far_h2.atoms.end());

  const BasisLibrary library = BasisLibrary::load("def2-SVP");
  const ScfResult a = run_rhf(hf, BasisSet::build(hf, library), tight());
  const ScfResult b = run_rhf(far_h2, BasisSet::build(far_h2, library), tight());
  const ScfResult ab = run_rhf(both, BasisSet::build(both, library), tight());
  REQUIRE(a.converged);
  REQUIRE(b.converged);
  REQUIRE(ab.converged);
  // Neutral molecules at 200 Bohr: only the dipole of HF interacts with the
  // quadrupole of H2, far below the tolerance.
  CHECK_THAT(ab.energy_total, WithinAbs(a.energy_total + b.energy_total, 1e-7));
  for (std::size_t k = 0; k < 3; ++k) CHECK_THAT(ab.dipole[k], WithinAbs(a.dipole[k] + b.dipole[k], 1e-6));
}

TEST_CASE("RHF: closed shells only", "[scf][rhf]") {
  Molecule hydrogen;
  hydrogen.atoms = {{1, {0.0, 0.0, 0.0}}};
  const BasisLibrary library = BasisLibrary::load("STO-3G");
  CHECK_THROWS_AS(run_rhf(hydrogen, BasisSet::build(hydrogen, library)), std::invalid_argument);

  // H- has two electrons but STO-3G offers one orbital: that is allowed (one
  // doubly occupied orbital). He with charge -2 would need two orbitals.
  hydrogen.charge = -1;
  CHECK(run_rhf(hydrogen, BasisSet::build(hydrogen, library)).converged);
  Molecule helium;
  helium.atoms = {{2, {0.0, 0.0, 0.0}}};
  helium.charge = -2;
  CHECK_THROWS_AS(run_rhf(helium, BasisSet::build(helium, library)), std::invalid_argument);

  // Water cation: odd number of electrons.
  Molecule water = load("h2o");
  water.charge = 1;
  CHECK_THROWS_AS(run_rhf(water, BasisSet::build(water, library)), std::invalid_argument);
}

TEST_CASE("RHF: not converged within the iteration limit is reported", "[scf][rhf]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("STO-3G"));
  ScfOptions options;
  options.max_iterations = 2;
  const ScfResult scf = run_rhf(water, basis, options);
  CHECK_FALSE(scf.converged);
  CHECK(scf.iterations == 2);
}

namespace {

// Wraps the direct builder and counts calls: the SCF must obtain J and K
// through the interface only.
class CountingJK final : public JKBuilder {
 public:
  explicit CountingJK(const BasisSet& basis) : inner_(basis) {}
  void build(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
             std::vector<Eigen::MatrixXd>* K) override {
    ++calls;
    inner_.build(densities, J, K);
  }
  void build_long_range_exchange(const std::vector<Eigen::MatrixXd>& densities, double omega,
                                 std::vector<Eigen::MatrixXd>& K) override {
    inner_.build_long_range_exchange(densities, omega, K);
  }
  int calls = 0;

 private:
  DirectJK inner_;
};

}  // namespace

TEST_CASE("RHF uses the J/K interface", "[scf][rhf]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("STO-3G"));
  CountingJK jk(basis);
  const ScfResult scf = run_rhf(water, basis, jk);
  CHECK(scf.converged);
  CHECK(jk.calls == scf.iterations);
  CHECK_THAT(scf.energy_total, WithinAbs(run_rhf(water, basis).energy_total, 1e-12));
}

TEST_CASE("DIIS extrapolation", "[scf][diis]") {
  const Eigen::MatrixXd F1 = Eigen::MatrixXd::Constant(2, 2, 1.0);
  const Eigen::MatrixXd F2 = Eigen::MatrixXd::Constant(2, 2, 4.0);
  Eigen::MatrixXd v(2, 2);
  v << 0.3, -0.2, 0.5, 0.1;

  SECTION("a single vector is returned unchanged") {
    Diis diis;
    CHECK(diis.extrapolate(F1, v) == F1);
    CHECK(diis.size() == 1);
  }
  SECTION("two antiparallel errors cancel exactly") {
    // e1 = v, e2 = -v/2: c1 + c2 = 1 and c1 - c2/2 = 0 give c1 = 1/3, c2 = 2/3.
    Diis diis;
    diis.extrapolate(F1, v);
    const Eigen::MatrixXd F = diis.extrapolate(F2, -0.5 * v);
    CHECK((F - (F1 / 3.0 + 2.0 * F2 / 3.0)).cwiseAbs().maxCoeff() < 1e-13);
  }
  SECTION("orthogonal errors of equal norm are averaged") {
    Eigen::MatrixXd e1 = Eigen::MatrixXd::Zero(2, 2), e2 = Eigen::MatrixXd::Zero(2, 2);
    e1(0, 0) = 1e-6;
    e2(1, 1) = 1e-6;
    Diis diis;
    diis.extrapolate(F1, e1);
    const Eigen::MatrixXd F = diis.extrapolate(F2, e2);
    CHECK((F - 0.5 * (F1 + F2)).cwiseAbs().maxCoeff() < 1e-12);
  }
  SECTION("the oldest vectors are dropped") {
    Diis diis(2);
    Eigen::MatrixXd e = v;
    for (int i = 0; i < 5; ++i) {
      e(0, 0) += 0.1 * (i + 1);
      diis.extrapolate(F1 * (i + 1), e);
      CHECK(diis.size() <= 2);
    }
  }
  SECTION("identical errors make the system singular; the result stays finite") {
    Diis diis;
    diis.extrapolate(F1, v);
    const Eigen::MatrixXd F = diis.extrapolate(F2, v);
    CHECK(F.allFinite());
  }
  CHECK_THROWS_AS(Diis(0), std::invalid_argument);
}
