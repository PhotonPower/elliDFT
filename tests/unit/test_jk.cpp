// Coulomb and exchange matrices from the direct builder against the explicit
// contraction of the full integral tensor.

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <vector>

#include "basis/basis_set.hpp"
#include "core/parallel.hpp"
#include "full_eri.hpp"
#include "scf/jk.hpp"

using namespace ellidft;

namespace {

const char* const test_basis = R"(
BASIS "ao basis" SPHERICAL
H    S
      2.9   0.3
      0.5   0.8
H    P
      0.8   1.0
O    S
      7.0   0.4
      1.1   0.7
O    S
      0.35  1.0
O    P
      1.9   0.5
      0.45  0.6
O    D
      1.0   1.0
END
)";

BasisSet test_basis_set(bool spherical) {
  Molecule mol;
  mol.atoms = {{8, {0.1, -0.2, 0.3}}, {1, {1.4, 0.9, -0.5}}, {1, {-0.9, 1.2, 1.1}}};
  std::istringstream in(test_basis);
  return BasisSet::build(mol, BasisLibrary::parse_nwchem(in, "test"), spherical);
}

// Deterministic pseudo-random matrix without symmetry.
Eigen::MatrixXd general_matrix(int n, unsigned seed) {
  Eigen::MatrixXd m(n, n);
  unsigned state = seed;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      state = state * 1664525u + 1013904223u;
      m(i, j) = static_cast<double>(state >> 8) / static_cast<double>(1u << 24) - 0.5;
    }
  return m;
}

// Definitions, term by term.
Eigen::MatrixXd coulomb_explicit(const testing::FullEri& eri, const Eigen::MatrixXd& D) {
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(eri.n, eri.n);
  for (int m = 0; m < eri.n; ++m)
    for (int n = 0; n < eri.n; ++n)
      for (int l = 0; l < eri.n; ++l)
        for (int s = 0; s < eri.n; ++s) J(m, n) += eri(m, n, l, s) * D(l, s);
  return J;
}

Eigen::MatrixXd exchange_explicit(const testing::FullEri& eri, const Eigen::MatrixXd& D) {
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(eri.n, eri.n);
  for (int m = 0; m < eri.n; ++m)
    for (int n = 0; n < eri.n; ++n)
      for (int l = 0; l < eri.n; ++l)
        for (int s = 0; s < eri.n; ++s) K(m, n) += eri(m, l, n, s) * D(l, s);
  return K;
}

}  // namespace

TEST_CASE("direct J and K equal the contraction of the full tensor", "[scf][jk]") {
  for (const bool spherical : {true, false}) {
    INFO((spherical ? "spherical" : "cartesian"));
    const BasisSet basis = test_basis_set(spherical);
    const int n = basis.n_functions();
    const testing::FullEri eri = testing::full_eri_tensor(basis);
    DirectJK jk(basis, DirectJKOptions{0.0});  // no screening: exact comparison

    const Eigen::MatrixXd general = general_matrix(n, 7);
    const Eigen::MatrixXd symmetric = general + general.transpose();
    const Eigen::MatrixXd antisymmetric = general - general.transpose();
    const std::vector<Eigen::MatrixXd> densities = {general, symmetric, antisymmetric};

    std::vector<Eigen::MatrixXd> J, K;
    jk.build(densities, &J, &K);
    REQUIRE(J.size() == 3);
    REQUIRE(K.size() == 3);
    for (std::size_t d = 0; d < densities.size(); ++d) {
      INFO("density " << d);
      CHECK((J[d] - coulomb_explicit(eri, densities[d])).cwiseAbs().maxCoeff() < 1e-12);
      CHECK((K[d] - exchange_explicit(eri, densities[d])).cwiseAbs().maxCoeff() < 1e-12);
    }
    // J sees only the symmetric part of the density; K of an antisymmetric
    // density is antisymmetric.
    CHECK(J[2].cwiseAbs().maxCoeff() < 1e-12);
    CHECK((K[2] + K[2].transpose()).cwiseAbs().maxCoeff() < 1e-12);
    CHECK((K[1] - K[1].transpose()).cwiseAbs().maxCoeff() < 1e-12);

    SECTION("only J or only K") {
      std::vector<Eigen::MatrixXd> J_only, K_only;
      jk.build({general}, &J_only, nullptr);
      jk.build({general}, nullptr, &K_only);
      CHECK((J_only[0] - J[0]).cwiseAbs().maxCoeff() < 1e-13);
      CHECK((K_only[0] - K[0]).cwiseAbs().maxCoeff() < 1e-13);
    }
  }
}

// Schwarz screening: with the default threshold the matrices change by far
// less than the threshold times the number of functions, and for a stretched
// system quartets are actually skipped.
TEST_CASE("direct J and K with Schwarz screening", "[scf][jk]") {
  Molecule mol;
  mol.atoms = {{8, {0.1, -0.2, 0.3}}, {1, {1.4, 0.9, -0.5}}, {1, {-0.9, 1.2, 1.1}}, {8, {14.0, 0.5, -0.4}},
               {1, {15.2, 1.6, 0.3}}, {1, {13.1, 1.9, -1.2}}};
  std::istringstream in(test_basis);
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::parse_nwchem(in, "test"), true);
  const Eigen::MatrixXd D = general_matrix(basis.n_functions(), 3);

  DirectJK exact(basis, DirectJKOptions{0.0}), screened(basis);
  std::vector<Eigen::MatrixXd> J0, K0, J1, K1;
  exact.build({D}, &J0, &K0);
  screened.build({D}, &J1, &K1);
  CHECK(exact.quartets_skipped() == 0);
  CHECK(screened.quartets_skipped() > 0);
  CHECK(screened.quartets_computed() + screened.quartets_skipped() == exact.quartets_computed());
  CHECK((J1[0] - J0[0]).cwiseAbs().maxCoeff() < 1e-10);
  CHECK((K1[0] - K0[0]).cwiseAbs().maxCoeff() < 1e-10);

  // A coarse threshold skips at least as much and stays within the bound
  // threshold * (number of function pairs).
  DirectJK coarse(basis, DirectJKOptions{1e-6});
  std::vector<Eigen::MatrixXd> J2, K2;
  coarse.build({D}, &J2, &K2);
  CHECK(coarse.quartets_skipped() >= screened.quartets_skipped());
  const double bound = 1e-6 * basis.n_functions() * basis.n_functions();
  CHECK((J2[0] - J0[0]).cwiseAbs().maxCoeff() < bound);
  CHECK((K2[0] - K0[0]).cwiseAbs().maxCoeff() < bound);
}

// Rule 11 of the specification: results depend on the number of threads only
// through the order of summation.
TEST_CASE("direct J and K do not depend on the number of threads", "[scf][jk][threads]") {
  const BasisSet basis = test_basis_set(true);
  const Eigen::MatrixXd D = general_matrix(basis.n_functions(), 11);
  const int available = max_threads();
  std::vector<Eigen::MatrixXd> J1, K1, Jn, Kn;
  set_num_threads(1);
  DirectJK(basis).build({D}, &J1, &K1);
  set_num_threads(std::max(available, 4));
  DirectJK(basis).build({D}, &Jn, &Kn);
  set_num_threads(available);
  CHECK((Jn[0] - J1[0]).cwiseAbs().maxCoeff() < 1e-12);
  CHECK((Kn[0] - K1[0]).cwiseAbs().maxCoeff() < 1e-12);
}

TEST_CASE("direct J and K reject densities of the wrong size", "[scf][jk]") {
  const BasisSet basis = test_basis_set(true);
  DirectJK jk(basis);
  std::vector<Eigen::MatrixXd> J, K;
  CHECK_THROWS_AS(jk.build({Eigen::MatrixXd::Zero(3, 3)}, &J, &K), std::invalid_argument);
}
