// Density fitting: checks that need no reference code.

#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <vector>

#include "core/xyz.hpp"
#include "df/df_integrals.hpp"
#include "df/df_jk.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;

namespace {

Molecule water() { return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2o.xyz"); }

Eigen::MatrixXd general_matrix(Eigen::Index n, unsigned seed) {
  Eigen::MatrixXd m(n, n);
  unsigned state = seed;
  for (Eigen::Index i = 0; i < n; ++i)
    for (Eigen::Index j = 0; j < n; ++j) {
      state = state * 1664525u + 1013904223u;
      m(i, j) = static_cast<double>(state >> 8) / static_cast<double>(1u << 24) - 0.5;
    }
  return m;
}

}  // namespace

// J and K from DensityFittedJK against the defining formula written out with
// dense linear algebra: (mu nu|lambda sigma) ~ sum_PQ (mu nu|P) [V^-1]_PQ (Q|lambda sigma).
TEST_CASE("density-fitted J and K equal the defining formula", "[df][jk]") {
  const Molecule mol = water();
  for (const char* basis_name : {"STO-3G", "6-31G*"}) {  // spherical and cartesian orbital basis
    INFO(basis_name);
    const BasisSet basis = BasisSet::build(mol, BasisLibrary::load(basis_name));
    const BasisSet aux = BasisSet::build(mol, BasisLibrary::load("def2-universal-JKFIT"));
    const Eigen::Index n = basis.n_functions();

    const Eigen::MatrixXd V = two_center_coulomb(aux);
    const Eigen::MatrixXd A = three_center_coulomb(basis, aux);
    const Eigen::MatrixXd Vinv = V.inverse();
    auto pair = [](Eigen::Index i, Eigen::Index j) { return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; };
    const Eigen::MatrixXd fitted_pairs = A * Vinv * A.transpose();  // (pair | pair)
    auto fitted = [&](Eigen::Index m, Eigen::Index nu, Eigen::Index l, Eigen::Index s) -> double {
      return fitted_pairs(pair(m, nu), pair(l, s));
    };

    const Eigen::MatrixXd D = general_matrix(n, 5);
    Eigen::MatrixXd J_ref = Eigen::MatrixXd::Zero(n, n), K_ref = Eigen::MatrixXd::Zero(n, n);
    for (Eigen::Index m = 0; m < n; ++m)
      for (Eigen::Index nu = 0; nu < n; ++nu)
        for (Eigen::Index l = 0; l < n; ++l)
          for (Eigen::Index s = 0; s < n; ++s) {
            J_ref(m, nu) += fitted(m, nu, l, s) * D(l, s);
            K_ref(m, nu) += fitted(m, l, nu, s) * D(l, s);
          }

    DensityFittedJK jk(basis, aux);
    CHECK(jk.n_auxiliary() == aux.n_functions());
    std::vector<Eigen::MatrixXd> J, K;
    jk.build({D, D.transpose()}, &J, &K);
    // The comparison values use the explicit inverse of the fitting metric,
    // which is ill conditioned; their rounding error is of the order 1e-9 and
    // depends on the linear algebra library (1.4e-9 seen with an external
    // BLAS). Hence 1e-8 here.
    CHECK((J[0] - J_ref).cwiseAbs().maxCoeff() < 1e-8);
    CHECK((K[0] - K_ref).cwiseAbs().maxCoeff() < 1e-8);
    CHECK((J[1] - J_ref).cwiseAbs().maxCoeff() < 1e-8);              // J sees the symmetric part
    CHECK((K[1] - K_ref.transpose()).cwiseAbs().maxCoeff() < 1e-8);  // K[D^T] = K[D]^T

    // The fitted matrices are close to the exact ones.
    DirectJK exact(basis, DirectJKOptions{0.0});
    std::vector<Eigen::MatrixXd> J0, K0;
    exact.build({D}, &J0, &K0);
    CHECK((J[0] - J0[0]).cwiseAbs().maxCoeff() < 2e-2);
    CHECK((K[0] - K0[0]).cwiseAbs().maxCoeff() < 2e-2);

    std::vector<Eigen::MatrixXd> J_only;
    jk.build({D}, &J_only, nullptr);
    CHECK((J_only[0] - J[0]).cwiseAbs().maxCoeff() < 1e-13);
    CHECK_THROWS_AS(jk.build({Eigen::MatrixXd::Zero(2, 2)}, &J, &K), std::invalid_argument);
  }
}

TEST_CASE("two-centre metric is symmetric and positive definite", "[df]") {
  const Molecule mol = water();
  const BasisSet aux = BasisSet::build(mol, BasisLibrary::load("def2-universal-JFIT"));
  const Eigen::MatrixXd V = two_center_coulomb(aux);
  CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-13);
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(V);
  CHECK(eigen.eigenvalues().minCoeff() > 0.0);
}

TEST_CASE("Hartree-Fock energy with density fitting is close to the exact one", "[df][scf]") {
  const Molecule mol = water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("6-31G*"));
  const BasisSet aux = BasisSet::build(mol, BasisLibrary::load("def2-universal-JKFIT"));
  DensityFittedJK jk(basis, aux);
  const ScfResult fitted = run_rhf(mol, basis, jk);
  const ScfResult exact = run_rhf(mol, basis);
  REQUIRE(fitted.converged);
  REQUIRE(exact.converged);
  const double error = fitted.energy_total - exact.energy_total;
  INFO("fitting error " << error);
  CHECK(std::abs(error) < 5e-4);
  CHECK(std::abs(error) > 1e-9);  // it is an approximation, not the same calculation
}
