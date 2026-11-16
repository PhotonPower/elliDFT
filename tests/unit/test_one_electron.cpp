// One-electron integrals: checks that need no reference code.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <sstream>
#include <vector>

#include "basis/basis_set.hpp"
#include "integrals/one_electron.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Artificial basis with all angular momenta up to g on two elements.
const char* const test_basis = R"(
BASIS "ao basis" SPHERICAL
H    S
      3.1   0.25
      0.6   0.80
H    P
      0.9   1.0
H    D
      0.8   1.0
O    S
      9.0   0.3
      1.8   0.7
O    S
      0.4   1.0
O    P
      2.2   0.4
      0.5   0.7
O    D
      1.1   1.0
O    F
      0.9   1.0
O    G
      0.7   1.0
END
)";

Molecule test_molecule(const Vec3& shift = {0.0, 0.0, 0.0}) {
  Molecule mol;
  mol.atoms = {{8, {0.1, -0.2, 0.3}}, {1, {1.4, 0.9, -0.5}}, {1, {-0.9, 1.2, 1.1}}};
  for (Atom& atom : mol.atoms)
    for (std::size_t k = 0; k < 3; ++k) atom.r[k] += shift[k];
  return mol;
}

BasisSet test_basis_set(const Molecule& mol, bool spherical) {
  std::istringstream in(test_basis);
  return BasisSet::build(mol, BasisLibrary::parse_nwchem(in, "test"), spherical);
}

double asymmetry(const Eigen::MatrixXd& m) { return (m - m.transpose()).cwiseAbs().maxCoeff(); }

}  // namespace

TEST_CASE("one-electron matrices are symmetric", "[integrals][one_electron]") {
  for (const bool spherical : {true, false}) {
    INFO((spherical ? "spherical" : "cartesian"));
    const Molecule mol = test_molecule();
    const BasisSet basis = test_basis_set(mol, spherical);
    // Both orders of every shell pair are computed independently, so this is
    // a genuine check.
    CHECK(asymmetry(overlap_matrix(basis)) < 1e-13);
    CHECK(asymmetry(kinetic_matrix(basis)) < 1e-12);
    CHECK(asymmetry(nuclear_attraction_matrix(basis, mol)) < 1e-12);
    for (const auto& d : dipole_matrices(basis)) CHECK(asymmetry(d) < 1e-12);
  }
}

TEST_CASE("overlap and kinetic matrices: definiteness and normalisation", "[integrals][one_electron]") {
  const Molecule mol = test_molecule();
  SECTION("spherical functions are normalised to 1") {
    const BasisSet basis = test_basis_set(mol, true);
    const Eigen::MatrixXd S = overlap_matrix(basis);
    CHECK((S.diagonal().array() - 1.0).abs().maxCoeff() < 1e-13);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(S), eig_T(kinetic_matrix(basis));
    CHECK(eig.eigenvalues().minCoeff() > 0.0);
    CHECK(eig_T.eigenvalues().minCoeff() > 0.0);
    // Functions of different angular momentum on the same centre are orthogonal.
    // Oxygen shells: s s p d f g -> functions 0, 1, 2-4, 5-9, 10-16, 17-25.
    CHECK(S.block(0, 2, 2, 24).cwiseAbs().maxCoeff() < 1e-13);    // s with p, d, f, g
    CHECK(S.block(2, 5, 3, 21).cwiseAbs().maxCoeff() < 1e-13);    // p with d, f, g
    CHECK(S.block(5, 10, 5, 16).cwiseAbs().maxCoeff() < 1e-13);   // d with f, g
    CHECK(S.block(10, 17, 7, 9).cwiseAbs().maxCoeff() < 1e-13);   // f with g
    // Within a shell the functions are orthonormal.
    CHECK((S.block(17, 17, 9, 9) - Eigen::MatrixXd::Identity(9, 9)).cwiseAbs().maxCoeff() < 1e-13);
  }
  SECTION("cartesian functions carry the radial normalisation only") {
    const BasisSet basis = test_basis_set(mol, false);
    const Eigen::MatrixXd S = overlap_matrix(basis);
    // Oxygen: s s p (0..4), d starts at 5: xx xy xz yy yz zz.
    // <xx|xx> = integral x^4 over the unit sphere = 4 pi / 5, <xy|xy> = 4 pi / 15.
    CHECK_THAT(S(5, 5), WithinRel(4.0 * M_PI / 5.0, 1e-13));
    CHECK_THAT(S(6, 6), WithinRel(4.0 * M_PI / 15.0, 1e-13));
    CHECK_THAT(S(5, 8), WithinRel(4.0 * M_PI / 15.0, 1e-13));  // <xx|yy>
    // s and p are normalised to 1.
    for (int i = 0; i < 5; ++i) CHECK_THAT(S(i, i), WithinRel(1.0, 1e-13));
  }
}

TEST_CASE("one-electron matrices under translation", "[integrals][one_electron]") {
  const Vec3 shift = {0.7, -1.3, 2.1};
  const Molecule mol = test_molecule(), moved = test_molecule(shift);
  for (const bool spherical : {true, false}) {
    INFO((spherical ? "spherical" : "cartesian"));
    const BasisSet basis = test_basis_set(mol, spherical), basis_moved = test_basis_set(moved, spherical);
    const Eigen::MatrixXd S = overlap_matrix(basis);
    CHECK((overlap_matrix(basis_moved) - S).cwiseAbs().maxCoeff() < 1e-13);
    CHECK((kinetic_matrix(basis_moved) - kinetic_matrix(basis)).cwiseAbs().maxCoeff() < 1e-12);
    CHECK((nuclear_attraction_matrix(basis_moved, moved) - nuclear_attraction_matrix(basis, mol))
              .cwiseAbs().maxCoeff() < 1e-11);
    // The dipole operator moves with the molecule if the origin moves too ...
    const auto d = dipole_matrices(basis);
    const auto d_moved_origin = dipole_matrices(basis_moved, shift);
    // ... and changes by shift * S if the origin stays.
    const auto d_moved = dipole_matrices(basis_moved);
    for (std::size_t k = 0; k < 3; ++k) {
      CHECK((d_moved_origin[k] - d[k]).cwiseAbs().maxCoeff() < 1e-12);
      CHECK((d_moved[k] - d[k] - shift[k] * S).cwiseAbs().maxCoeff() < 1e-12);
    }
  }
}

TEST_CASE("overlap through the Hermite expansion of shell pairs", "[integrals][one_electron]") {
  for (const bool spherical : {true, false}) {
    const BasisSet basis = test_basis_set(test_molecule(), spherical);
    CHECK((overlap_matrix_from_shell_pairs(basis) - overlap_matrix(basis)).cwiseAbs().maxCoeff() < 1e-13);
  }
}

// Closed forms for s functions, from the Gaussian product rule: two
// normalised s Gaussians with exponents a, b at distance R have
//   S = (2 sqrt(ab)/(a+b))^(3/2) exp(-ab/(a+b) R^2),
//   T = ab/(a+b) (3 - 2 ab/(a+b) R^2) S.
TEST_CASE("s functions: closed forms for overlap and kinetic energy", "[integrals][one_electron]") {
  const double a = 0.9, b = 1.6, R = 1.3;
  std::vector<Shell> shells(2);
  shells[0].atom = 0;
  shells[0].center = {0.0, 0.0, 0.0};
  shells[0].exponents = {a};
  shells[0].coefficients = {1.0};
  shells[1].atom = 1;
  shells[1].center = {0.0, 0.0, R};
  shells[1].exponents = {b};
  shells[1].coefficients = {1.0};
  const BasisSet basis = BasisSet::from_shells("s", true, shells);
  const double mu = a * b / (a + b);
  const double S01 = std::pow(2.0 * std::sqrt(a * b) / (a + b), 1.5) * std::exp(-mu * R * R);
  const Eigen::MatrixXd S = overlap_matrix(basis), T = kinetic_matrix(basis);
  CHECK_THAT(S(0, 0), WithinRel(1.0, 1e-14));
  CHECK_THAT(S(0, 1), WithinRel(S01, 1e-14));
  CHECK_THAT(T(0, 0), WithinRel(1.5 * a, 1e-14));
  CHECK_THAT(T(0, 1), WithinRel(mu * (3.0 - 2.0 * mu * R * R) * S01, 1e-13));

  // Nuclear attraction of a normalised s Gaussian to a unit charge at its own
  // centre: -2 sqrt(2a/pi).
  Molecule mol;
  mol.atoms = {{1, {0.0, 0.0, 0.0}}};
  const BasisSet one = BasisSet::from_shells("s", true, {shells[0]});
  CHECK_THAT(nuclear_attraction_matrix(one, mol)(0, 0), WithinRel(-2.0 * std::sqrt(2.0 * a / M_PI), 1e-14));

  // Dipole: the product of the two s functions is centred at P = b R / (a+b).
  const auto d = dipole_matrices(basis);
  CHECK_THAT(d[2](0, 1), WithinRel(b * R / (a + b) * S01, 1e-13));
  CHECK_THAT(d[0](0, 1), WithinAbs(0.0, 1e-15));
}
