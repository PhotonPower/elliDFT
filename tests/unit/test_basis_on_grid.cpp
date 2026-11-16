// Basis functions on a grid: checks that need no reference code.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <sstream>
#include <vector>

#include "core/xyz.hpp"
#include "grid/basis_on_grid.hpp"
#include "grid/molecular_grid.hpp"
#include "integrals/one_electron.hpp"

using namespace ellidft;

namespace {

const char* const test_basis = R"(
BASIS "ao basis" SPHERICAL
H    S
      3.1   0.25
      0.6   0.80
H    P
      0.9   1.0
O    S
      9.0   0.3
      1.8   0.7
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

Molecule test_molecule() {
  Molecule mol;
  mol.atoms = {{8, {0.1, -0.2, 0.3}}, {1, {1.4, 0.9, -0.5}}, {1, {-0.9, 1.2, 1.1}}};
  return mol;
}

BasisSet test_basis_set(bool spherical) {
  std::istringstream in(test_basis);
  return BasisSet::build(test_molecule(), BasisLibrary::parse_nwchem(in, "test"), spherical);
}

}  // namespace

TEST_CASE("gradients of basis functions agree with finite differences", "[grid][basis_on_grid]") {
  const std::vector<Vec3> points = {{0.2, 0.1, -0.4}, {1.3, 0.8, -0.3}, {-2.0, 1.5, 0.7}, {0.1, -0.2, 0.3}};
  const double h = 1e-4;
  for (const bool spherical : {true, false}) {
    INFO((spherical ? "spherical" : "cartesian"));
    const BasisSet basis = test_basis_set(spherical);
    BasisValues center;
    evaluate_basis(basis, points.data(), points.size(), true, center, 0.0);
    REQUIRE(center.value.rows() == 4);
    REQUIRE(center.value.cols() == basis.n_functions());

    const Eigen::MatrixXd* analytic[3] = {&center.dx, &center.dy, &center.dz};
    for (std::size_t axis = 0; axis < 3; ++axis) {
      // Central differences at two step sizes with Richardson extrapolation.
      auto difference = [&](double step) {
        std::vector<Vec3> plus = points, minus = points;
        for (std::size_t p = 0; p < points.size(); ++p) {
          plus[p][axis] += step;
          minus[p][axis] -= step;
        }
        BasisValues a, b;
        evaluate_basis(basis, plus.data(), plus.size(), false, a, 0.0);
        evaluate_basis(basis, minus.data(), minus.size(), false, b, 0.0);
        return Eigen::MatrixXd((a.value - b.value) / (2.0 * step));
      };
      const Eigen::MatrixXd numeric = (4.0 * difference(h) - difference(2.0 * h)) / 3.0;
      CHECK((numeric - *analytic[axis]).cwiseAbs().maxCoeff() < 1e-9);
    }
  }
}

TEST_CASE("values only: gradients are not allocated", "[grid][basis_on_grid]") {
  const BasisSet basis = test_basis_set(true);
  const std::vector<Vec3> points = {{0.2, 0.1, -0.4}};
  BasisValues values;
  evaluate_basis(basis, points.data(), points.size(), false, values);
  CHECK(values.value.rows() == 1);
  CHECK(values.dx.size() == 0);
}

// Ties grid and basis functions together: the overlap matrix by quadrature
// equals the analytic one.
TEST_CASE("overlap matrix by quadrature on the molecular grid", "[grid][basis_on_grid]") {
  const Molecule water = read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2o.xyz");
  for (const char* basis_name : {"def2-SVP", "6-31G*"}) {
    INFO(basis_name);
    const BasisSet basis = BasisSet::build(water, BasisLibrary::load(basis_name));
    const MolecularGrid grid = build_molecular_grid(water, {75, 302});
    const Eigen::Index n = basis.n_functions();
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
    BasisValues values;
    const std::size_t block = 512;
    for (std::size_t start = 0; start < grid.size(); start += block) {
      const std::size_t count = std::min(block, grid.size() - start);
      evaluate_basis(basis, &grid.points[start], count, false, values);
      for (std::size_t p = 0; p < count; ++p) {
        const auto row = values.value.row(static_cast<Eigen::Index>(p));
        S.noalias() += grid.weights[start + p] * row.transpose() * row;
      }
    }
    CHECK((S - overlap_matrix(basis)).cwiseAbs().maxCoeff() < 2e-6);
  }
}

TEST_CASE("distant shells are skipped only where they are negligible", "[grid][basis_on_grid]") {
  const BasisSet basis = test_basis_set(true);
  // A block of points 40 Bohr away: every function is below 1e-16 there.
  std::vector<Vec3> far;
  for (int i = 0; i < 5; ++i) far.push_back({40.0 + 0.1 * i, 0.0, 0.0});
  BasisValues screened, exact;
  evaluate_basis(basis, far.data(), far.size(), true, screened);
  evaluate_basis(basis, far.data(), far.size(), true, exact, 0.0);
  CHECK(screened.value.cwiseAbs().maxCoeff() == 0.0);
  CHECK(exact.value.cwiseAbs().maxCoeff() < 1e-16);
  // Near the molecule nothing is lost.
  const std::vector<Vec3> near = {{0.5, 0.5, 0.5}, {3.0, -2.0, 1.0}, {6.0, 5.0, -4.0}};
  evaluate_basis(basis, near.data(), near.size(), true, screened);
  evaluate_basis(basis, near.data(), near.size(), true, exact, 0.0);
  CHECK((screened.value - exact.value).cwiseAbs().maxCoeff() < 1e-15);
  CHECK((screened.dz - exact.dz).cwiseAbs().maxCoeff() < 1e-15);
}
