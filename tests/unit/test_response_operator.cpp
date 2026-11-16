// Acceptance criterion M5: for small molecules the explicitly built matrices
// A and B agree with the matrix-free application to 1e-9.

#include <catch2/catch_test_macros.hpp>

#include <string>

#include "core/xyz.hpp"
#include "explicit_ab.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;

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

}  // namespace

TEST_CASE("explicit A and B equal the matrix-free application", "[response][m5]") {
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("STO-3G"));
  for (const char* name : {"HF", "LDA", "PBE", "PBE0", "CAM-B3LYP"}) {
    DYNAMIC_SECTION(name) {
      Method method;
      method.functional = Functional::by_name(name);
      method.grid = {30, 110};
      DirectJK jk(basis, DirectJKOptions{0.0});
      const ScfResult scf = run_scf(water, basis, method, jk, tight());
      REQUIRE(scf.converged);

      const testing::ExplicitAB ab = testing::explicit_ab(water, basis, method, scf);
      const ResponseOperator op(water, basis, method, scf, jk);
      REQUIRE(op.size() == ab.A.rows());
      REQUIRE(op.n_occ() == 5);
      REQUIRE(op.n_virt() == 2);

      Eigen::MatrixXd plus, minus;
      op.apply(Eigen::MatrixXd::Identity(op.size(), op.size()), &plus, &minus);
      const double dev_plus = (plus - (ab.A + ab.B)).cwiseAbs().maxCoeff();
      const double dev_minus = (minus - (ab.A - ab.B)).cwiseAbs().maxCoeff();
      INFO("A+B: " << dev_plus << ", A-B: " << dev_minus);
      CHECK(dev_plus <= 1e-9);
      CHECK(dev_minus <= 1e-9);
      SUCCEED(name << ": A+B " << dev_plus << ", A-B " << dev_minus);

      // Symmetry of A and B; A-B is diagonal without exact exchange.
      CHECK((ab.A - ab.A.transpose()).cwiseAbs().maxCoeff() < 1e-10);
      CHECK((ab.B - ab.B.transpose()).cwiseAbs().maxCoeff() < 1e-10);
      if (method.functional.exact_exchange == 0.0 && method.functional.long_range_exchange == 0.0) {
        const Eigen::MatrixXd diagonal = op.orbital_energy_differences().asDiagonal();
        CHECK((minus - diagonal).cwiseAbs().maxCoeff() < 1e-12);
      }
      // A stable ground state: A+B and A-B are positive definite.
      CHECK(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(ab.A + ab.B).eigenvalues().minCoeff() > 0.0);
      CHECK(Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>(ab.A - ab.B).eigenvalues().minCoeff() > 0.0);

      // Only one of the two products requested.
      Eigen::MatrixXd only_plus, only_minus;
      const Eigen::MatrixXd v = Eigen::MatrixXd::Identity(op.size(), 2);
      op.apply(v, &only_plus, nullptr);
      op.apply(v, nullptr, &only_minus);
      CHECK((only_plus - plus.leftCols(2)).cwiseAbs().maxCoeff() < 1e-12);
      CHECK((only_minus - minus.leftCols(2)).cwiseAbs().maxCoeff() < 1e-12);
    }
  }
}

TEST_CASE("matrix-free application with density fitting equals explicit fitted matrices approximately",
          "[response]") {
  // With a larger basis (virtual space of 14 orbitals) and cartesian d
  // functions; exact integrals.
  const Molecule water = load("h2o");
  const BasisSet basis = BasisSet::build(water, BasisLibrary::load("6-31G*"));
  Method method;
  method.functional = Functional::by_name("PBE0");
  method.grid = {20, 50};
  DirectJK jk(basis, DirectJKOptions{0.0});
  const ScfResult scf = run_scf(water, basis, method, jk, tight());
  REQUIRE(scf.converged);
  const testing::ExplicitAB ab = testing::explicit_ab(water, basis, method, scf);
  const ResponseOperator op(water, basis, method, scf, jk);
  REQUIRE(op.size() == 5 * 14);
  Eigen::MatrixXd plus, minus;
  op.apply(Eigen::MatrixXd::Identity(op.size(), op.size()), &plus, &minus);
  CHECK((plus - (ab.A + ab.B)).cwiseAbs().maxCoeff() <= 1e-9);
  CHECK((minus - (ab.A - ab.B)).cwiseAbs().maxCoeff() <= 1e-9);
  CHECK_THROWS_AS(op.apply(Eigen::MatrixXd::Zero(3, 1), &plus, &minus), std::invalid_argument);
  CHECK(op.dipole_vectors().rows() == op.size());
  CHECK(op.dipole_vectors().cols() == 3);
}
