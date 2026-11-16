// Lebedev grids. Acceptance criterion M2: the grids integrate spherical
// harmonics up to their nominal degree to 1e-12.
//
// The test uses orthonormal real spherical harmonics. Monomials would not do:
// a monomial of degree L contains the harmonics of degree L only with a
// weight of the order 2^-L, so it cannot reveal a grid that fails there.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "grid/lebedev.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;

namespace {

// Orthonormal real spherical harmonics up to degree lmax at the unit vector
// (x, y, z): values[l*l + l + m], m < 0 sine type, m > 0 cosine type.
//
// Theta_lm(u) = N_lm P_l^m(u), u = cos(theta), normalised to
// integral_{-1}^{1} Theta_lm^2 du = 1, from the relations of the associated
// Legendre functions as given in Wikipedia, "Associated Legendre
// polynomials" (consulted 2026-11-14):
//   P_l^l = (-1)^l (2l-1)!! (1-u^2)^(l/2),  P_{l+1}^l = u (2l+1) P_l^l,
//   (l-m+1) P_{l+1}^m = (2l+1) u P_l^m - (l+m) P_{l-1}^m,
//   integral P_l^m P_l^m du = 2 (l+m)! / ((2l+1) (l-m)!).
// Dividing by the norms gives the recurrences below, free of factorials.
// Y_l0 = Theta_l0 / sqrt(2 pi), Y_lm = Theta_l|m| cos(m phi) / sqrt(pi) or
// sin(|m| phi) / sqrt(pi). The normalisation is checked in the test itself.
void real_spherical_harmonics(int lmax, double x, double y, double z, std::vector<double>& values) {
  values.assign(static_cast<std::size_t>((lmax + 1) * (lmax + 1)), 0.0);
  const double s = std::sqrt(std::max(0.0, 1.0 - z * z));
  const double phi = std::atan2(y, x);
  std::vector<double> theta(static_cast<std::size_t>(lmax + 1));
  double theta_mm = std::sqrt(0.5);  // Theta_00
  for (int m = 0; m <= lmax; ++m) {
    if (m > 0) theta_mm *= -s * std::sqrt((2.0 * m + 1.0) / (2.0 * m));
    theta[static_cast<std::size_t>(m)] = theta_mm;
    if (m < lmax) theta[static_cast<std::size_t>(m + 1)] = z * std::sqrt(2.0 * m + 3.0) * theta_mm;
    for (int l = m + 1; l < lmax; ++l) {
      const double a = std::sqrt((2.0 * l + 3.0) * (2.0 * l + 1.0) / ((l + 1.0 - m) * (l + 1.0 + m)));
      const double b = std::sqrt((2.0 * l + 3.0) * (l - m) * (l + m) / ((2.0 * l - 1.0) * (l + 1.0 - m) * (l + 1.0 + m)));
      theta[static_cast<std::size_t>(l + 1)] =
          a * z * theta[static_cast<std::size_t>(l)] - b * theta[static_cast<std::size_t>(l - 1)];
    }
    for (int l = m; l <= lmax; ++l) {
      const double value = theta[static_cast<std::size_t>(l)];
      const auto center = static_cast<std::size_t>(l * l + l);
      if (m == 0) {
        values[center] = value / std::sqrt(2.0 * M_PI);
      } else {
        values[center + static_cast<std::size_t>(m)] = value * std::cos(m * phi) / std::sqrt(M_PI);
        values[center - static_cast<std::size_t>(m)] = value * std::sin(m * phi) / std::sqrt(M_PI);
      }
    }
  }
}

// errors[l] = largest |sum_i w_i Y_lm(p_i) - (exact mean of Y_lm)| over m.
// The mean over the sphere is 1/sqrt(4 pi) for l = 0 and 0 otherwise.
std::vector<double> harmonic_errors(const AngularGrid& grid, int lmax) {
  std::vector<double> sums(static_cast<std::size_t>((lmax + 1) * (lmax + 1)), 0.0), values;
  for (std::size_t i = 0; i < grid.points.size(); ++i) {
    real_spherical_harmonics(lmax, grid.points[i][0], grid.points[i][1], grid.points[i][2], values);
    for (std::size_t k = 0; k < sums.size(); ++k) sums[k] += grid.weights[i] * values[k];
  }
  sums[0] -= 1.0 / std::sqrt(4.0 * M_PI);
  std::vector<double> errors(static_cast<std::size_t>(lmax + 1), 0.0);
  for (int l = 0; l <= lmax; ++l)
    for (int k = l * l; k < (l + 1) * (l + 1); ++k)
      errors[static_cast<std::size_t>(l)] = std::max(errors[static_cast<std::size_t>(l)], std::abs(sums[static_cast<std::size_t>(k)]));
  return errors;
}

void check_grids(int min_points, int max_points) {
  for (const int n_points : lebedev_sizes()) {
    if (n_points < min_points || n_points > max_points) continue;
    const AngularGrid grid = lebedev_grid(n_points);
    INFO(n_points << " points, degree " << grid.degree);
    REQUIRE(grid.points.size() == static_cast<std::size_t>(n_points));
    REQUIRE(grid.weights.size() == grid.points.size());
    CHECK(grid.degree == lebedev_degree(n_points));

    double weight_sum = 0.0, norm_error = 0.0;
    for (std::size_t i = 0; i < grid.points.size(); ++i) {
      weight_sum += grid.weights[i];
      const auto& p = grid.points[i];
      norm_error = std::max(norm_error, std::abs(p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - 1.0));
    }
    CHECK_THAT(weight_sum, WithinAbs(1.0, 1e-13));
    CHECK(norm_error < 1e-14);

    // Acceptance criterion: all harmonics up to the nominal degree.
    const std::vector<double> errors = harmonic_errors(grid, grid.degree + 1);
    double worst = 0.0;
    for (int l = 0; l <= grid.degree; ++l) worst = std::max(worst, errors[static_cast<std::size_t>(l)]);
    CHECK(worst <= 1e-12);
    // The nominal degree is sharp: one degree further the grid is not exact.
    CHECK(errors[static_cast<std::size_t>(grid.degree + 1)] > 1e-8);
    SUCCEED(n_points << " points, degree " << grid.degree << ": largest error " << worst << ", at degree + 1: "
                     << errors[static_cast<std::size_t>(grid.degree + 1)]);
  }
}

}  // namespace

// The harmonics used above are orthonormal: checked with a grid that
// integrates their products exactly (degree 41 covers l + l' <= 41).
TEST_CASE("test harmonics are orthonormal on the sphere", "[grid][lebedev]") {
  const AngularGrid grid = lebedev_grid(590);
  const int lmax = 20;
  const auto n = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
  std::vector<double> gram(n * n, 0.0), values;
  for (std::size_t i = 0; i < grid.points.size(); ++i) {
    real_spherical_harmonics(lmax, grid.points[i][0], grid.points[i][1], grid.points[i][2], values);
    for (std::size_t a = 0; a < n; ++a)
      for (std::size_t b = 0; b < n; ++b) gram[a * n + b] += 4.0 * M_PI * grid.weights[i] * values[a] * values[b];
  }
  double worst = 0.0;
  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t b = 0; b < n; ++b) worst = std::max(worst, std::abs(gram[a * n + b] - (a == b ? 1.0 : 0.0)));
  CHECK(worst < 1e-12);
  // First functions in closed form: Y_00 = 1/sqrt(4 pi), Y_10 = sqrt(3/(4 pi)) z.
  real_spherical_harmonics(2, 0.3, -0.4, std::sqrt(0.75), values);
  CHECK_THAT(values[0], WithinAbs(1.0 / std::sqrt(4.0 * M_PI), 1e-15));
  CHECK_THAT(values[2], WithinAbs(std::sqrt(3.0 / (4.0 * M_PI)) * std::sqrt(0.75), 1e-15));
}

// All 32 grids, 6 to 5810 points (degrees 3 to 131); about two seconds.
TEST_CASE("Lebedev grids integrate spherical harmonics up to the nominal degree", "[grid][lebedev][m2]") {
  check_grids(0, 5810);
}

TEST_CASE("Lebedev grids: list of sizes and invalid requests", "[grid][lebedev]") {
  const auto& sizes = lebedev_sizes();
  REQUIRE(sizes.size() == 32);
  CHECK(sizes.front() == 6);
  CHECK(sizes.back() == 5810);
  CHECK(std::is_sorted(sizes.begin(), sizes.end()));
  CHECK(lebedev_degree(302) == 29);
  CHECK(lebedev_degree(590) == 41);
  CHECK_THROWS_AS(lebedev_grid(100), std::invalid_argument);
  CHECK_THROWS_AS(lebedev_degree(0), std::invalid_argument);

  // Octahedral symmetry: with every point the grid contains its mirror
  // images and coordinate permutations, with the same weight.
  const AngularGrid grid = lebedev_grid(110);
  auto contains = [&](double x, double y, double z, double w) {
    for (std::size_t i = 0; i < grid.points.size(); ++i)
      if (grid.points[i][0] == x && grid.points[i][1] == y && grid.points[i][2] == z && grid.weights[i] == w)
        return true;
    return false;
  };
  for (std::size_t i = 0; i < grid.points.size(); ++i) {
    const auto& p = grid.points[i];
    const double w = grid.weights[i];
    CHECK(contains(-p[0] + 0.0, p[1], p[2], w));
    CHECK(contains(p[1], p[0], p[2], w));
    CHECK(contains(p[2], p[1], p[0], w));
  }
}
