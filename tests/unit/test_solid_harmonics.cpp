#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <map>
#include <tuple>
#include <vector>

#include "basis/basis_library.hpp"
#include "basis/basis_set.hpp"
#include "basis/solid_harmonics.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

using Powers = std::tuple<int, int, int>;
using Polynomial = std::map<Powers, double>;  // only non-zero coefficients

Polynomial as_polynomial(int l, const std::vector<double>& coefficients) {
  Polynomial p;
  const auto& powers = cartesian_powers(l);
  for (std::size_t c = 0; c < powers.size(); ++c)
    if (std::abs(coefficients[c]) > 1e-12) p[{powers[c].x, powers[c].y, powers[c].z}] = coefficients[c];
  return p;
}

void check_polynomial(int l, int m, const Polynomial& expected) {
  INFO("(l, m) = (" << l << ", " << m << ")");
  const Polynomial actual = as_polynomial(l, real_solid_harmonic(l, m));
  REQUIRE(actual.size() == expected.size());
  for (const auto& [powers, value] : expected) {
    INFO("monomial x^" << std::get<0>(powers) << " y^" << std::get<1>(powers) << " z^"
                       << std::get<2>(powers));
    REQUIRE(actual.count(powers) == 1);
    CHECK_THAT(actual.at(powers), WithinRel(value, 1e-13));
  }
}

// Laplacian of a homogeneous polynomial given by coefficients in
// cartesian_powers(l) order, as polynomial of degree l - 2.
Polynomial laplacian(int l, const std::vector<double>& coefficients) {
  Polynomial out;
  const auto& powers = cartesian_powers(l);
  for (std::size_t c = 0; c < powers.size(); ++c) {
    const auto [x, y, z] = std::tuple{powers[c].x, powers[c].y, powers[c].z};
    if (x >= 2) out[{x - 2, y, z}] += coefficients[c] * x * (x - 1);
    if (y >= 2) out[{x, y - 2, z}] += coefficients[c] * y * (y - 1);
    if (z >= 2) out[{x, y, z - 2}] += coefficients[c] * z * (z - 1);
  }
  return out;
}

// Numerical integral over the unit sphere, independent of the Gamma function
// formula. Variables u = cos(theta) and phi: composite Simpson rule in u (the
// integrand is smooth in u; error of order h^4) and the midpoint rule in phi
// (exact for trigonometric polynomials of degree below the number of points).
// A plain midpoint rule in theta is only accurate to about 1e-5 here.
double sphere_integral_numeric(int a, int b, int c) {
  const int n_u = 4000, n_phi = 64;
  double sum = 0.0;
  for (int i = 0; i <= n_u; ++i) {
    const double u = -1.0 + 2.0 * i / n_u;
    const double weight_u = (i == 0 || i == n_u) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - u * u));
    double ring = 0.0;
    for (int j = 0; j < n_phi; ++j) {
      const double phi = 2.0 * M_PI * (j + 0.5) / n_phi;
      ring += std::pow(sin_theta * std::cos(phi), a) * std::pow(sin_theta * std::sin(phi), b);
    }
    sum += weight_u * ring * std::pow(u, c);
  }
  return sum * (2.0 / n_u / 3.0) * (2.0 * M_PI / n_phi);
}

}  // namespace

TEST_CASE("cartesian monomials in lexicographic order", "[basis][cartesian]") {
  const auto& d = cartesian_powers(2);
  REQUIRE(d.size() == 6);
  // xx xy xz yy yz zz
  const int expected[6][3] = {{2, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, 2, 0}, {0, 1, 1}, {0, 0, 2}};
  for (std::size_t c = 0; c < 6; ++c) {
    CHECK(d[c].x == expected[c][0]);
    CHECK(d[c].y == expected[c][1]);
    CHECK(d[c].z == expected[c][2]);
  }
  for (int l = 0; l <= max_angular_momentum; ++l) {
    const auto& powers = cartesian_powers(l);
    REQUIRE(powers.size() == static_cast<std::size_t>(n_cartesian(l)));
    for (std::size_t c = 0; c < powers.size(); ++c) {
      CHECK(powers[c].x + powers[c].y + powers[c].z == l);
      CHECK(cartesian_index(powers[c].x, powers[c].y, powers[c].z) == static_cast<int>(c));
    }
  }
  CHECK_THROWS_AS(cartesian_powers(-1), std::out_of_range);
  CHECK_THROWS_AS(cartesian_powers(max_angular_momentum + 1), std::out_of_range);
}

TEST_CASE("integrals of monomials over the unit sphere", "[basis][solid_harmonics]") {
  CHECK_THAT(sphere_monomial_integral(0, 0, 0), WithinRel(4.0 * M_PI, 1e-15));
  CHECK(sphere_monomial_integral(1, 0, 0) == 0.0);
  CHECK(sphere_monomial_integral(2, 1, 0) == 0.0);
  CHECK(sphere_monomial_integral(0, 2, 3) == 0.0);
  for (const auto& [a, b, c] : std::vector<Powers>{{2, 0, 0}, {0, 0, 2}, {2, 2, 0}, {4, 0, 0},
                                                   {2, 2, 2}, {6, 4, 2}, {0, 12, 0}, {4, 4, 4}}) {
    INFO("x^" << a << " y^" << b << " z^" << c);
    CHECK_THAT(sphere_monomial_integral(a, b, c), WithinRel(sphere_integral_numeric(a, b, c), 1e-9));
  }
}

// Literature values: Ribaldone & Desmarais 2025 (arXiv:2412.16733v2),
// Tables I and II, coefficients D^{lm}_{tuv} of the unnormalised real solid
// harmonics; the triplets are the exponents (t, u, v) of x, y, z.
TEST_CASE("real solid harmonics equal the tabulated polynomials", "[basis][solid_harmonics]") {
  check_polynomial(0, 0, {{{0, 0, 0}, 1.0}});
  check_polynomial(1, -1, {{{0, 1, 0}, 1.0}});
  check_polynomial(1, 0, {{{0, 0, 1}, 1.0}});
  check_polynomial(1, 1, {{{1, 0, 0}, 1.0}});
  check_polynomial(2, -2, {{{1, 1, 0}, 6.0}});
  check_polynomial(2, -1, {{{0, 1, 1}, 3.0}});
  check_polynomial(2, 0, {{{0, 0, 2}, 1.0}, {{0, 2, 0}, -0.5}, {{2, 0, 0}, -0.5}});
  check_polynomial(2, 1, {{{1, 0, 1}, 3.0}});
  check_polynomial(2, 2, {{{0, 2, 0}, -3.0}, {{2, 0, 0}, 3.0}});
  check_polynomial(3, -3, {{{0, 3, 0}, -15.0}, {{2, 1, 0}, 45.0}});
  check_polynomial(3, -2, {{{1, 1, 1}, 30.0}});
  check_polynomial(3, -1, {{{0, 1, 2}, 6.0}, {{0, 3, 0}, -1.5}, {{2, 1, 0}, -1.5}});
  check_polynomial(3, 0, {{{0, 0, 3}, 1.0}, {{0, 2, 1}, -1.5}, {{2, 0, 1}, -1.5}});
  check_polynomial(3, 1, {{{1, 0, 2}, 6.0}, {{1, 2, 0}, -1.5}, {{3, 0, 0}, -1.5}});
  check_polynomial(3, 2, {{{0, 2, 1}, -15.0}, {{2, 0, 1}, 15.0}});
  check_polynomial(3, 3, {{{1, 2, 0}, -45.0}, {{3, 0, 0}, 15.0}});
  check_polynomial(4, -4, {{{1, 3, 0}, -420.0}, {{3, 1, 0}, 420.0}});
  check_polynomial(4, -3, {{{0, 3, 1}, -105.0}, {{2, 1, 1}, 315.0}});
  check_polynomial(4, -2, {{{1, 1, 2}, 90.0}, {{1, 3, 0}, -15.0}, {{3, 1, 0}, -15.0}});
  check_polynomial(4, -1, {{{0, 1, 3}, 10.0}, {{0, 3, 1}, -7.5}, {{2, 1, 1}, -7.5}});
  check_polynomial(4, 0, {{{0, 0, 4}, 1.0}, {{0, 2, 2}, -3.0}, {{0, 4, 0}, 0.375},
                          {{2, 0, 2}, -3.0}, {{2, 2, 0}, 0.75}, {{4, 0, 0}, 0.375}});
  check_polynomial(4, 1, {{{1, 0, 3}, 10.0}, {{1, 2, 1}, -7.5}, {{3, 0, 1}, -7.5}});
  check_polynomial(4, 2, {{{0, 2, 2}, -45.0}, {{0, 4, 0}, 7.5}, {{2, 0, 2}, 45.0}, {{4, 0, 0}, -7.5}});
  check_polynomial(4, 3, {{{1, 2, 1}, -315.0}, {{3, 0, 1}, 105.0}});
  check_polynomial(4, 4, {{{0, 4, 0}, 105.0}, {{2, 2, 0}, -630.0}, {{4, 0, 0}, 105.0}});
  // l = 5, Table II (selection covering both signs of m and m = 0)
  check_polynomial(5, -5, {{{0, 5, 0}, 945.0}, {{2, 3, 0}, -9450.0}, {{4, 1, 0}, 4725.0}});
  check_polynomial(5, -4, {{{1, 3, 1}, -3780.0}, {{3, 1, 1}, 3780.0}});
  check_polynomial(5, -2, {{{1, 1, 3}, 210.0}, {{1, 3, 1}, -105.0}, {{3, 1, 1}, -105.0}});
  check_polynomial(5, 0, {{{0, 0, 5}, 1.0}, {{0, 2, 3}, -5.0}, {{0, 4, 1}, 1.875},
                          {{2, 0, 3}, -5.0}, {{2, 2, 1}, 3.75}, {{4, 0, 1}, 1.875}});
  check_polynomial(5, 2, {{{0, 2, 3}, -105.0}, {{0, 4, 1}, 52.5}, {{2, 0, 3}, 105.0}, {{4, 0, 1}, -52.5}});
  check_polynomial(5, 4, {{{0, 4, 1}, 945.0}, {{2, 2, 1}, -5670.0}, {{4, 0, 1}, 945.0}});
  check_polynomial(5, 5, {{{1, 4, 0}, 4725.0}, {{3, 2, 0}, -9450.0}, {{5, 0, 0}, 945.0}});
  // l = 6, Table II (selection)
  check_polynomial(6, -6, {{{1, 5, 0}, 62370.0}, {{3, 3, 0}, -207900.0}, {{5, 1, 0}, 62370.0}});
  check_polynomial(6, -5, {{{0, 5, 1}, 10395.0}, {{2, 3, 1}, -103950.0}, {{4, 1, 1}, 51975.0}});
  check_polynomial(6, -3, {{{0, 3, 3}, -1260.0}, {{0, 5, 1}, 472.5}, {{2, 1, 3}, 3780.0},
                           {{2, 3, 1}, -945.0}, {{4, 1, 1}, -1417.5}});
  check_polynomial(6, 0, {{{0, 0, 6}, 1.0}, {{0, 2, 4}, -7.5}, {{0, 4, 2}, 5.625},
                          {{0, 6, 0}, -0.3125}, {{2, 0, 4}, -7.5}, {{2, 2, 2}, 11.25},
                          {{2, 4, 0}, -0.9375}, {{4, 0, 2}, 5.625}, {{4, 2, 0}, -0.9375},
                          {{6, 0, 0}, -0.3125}});
  check_polynomial(6, 3, {{{1, 2, 3}, -3780.0}, {{1, 4, 1}, 1417.5}, {{3, 0, 3}, 1260.0},
                          {{3, 2, 1}, 945.0}, {{5, 0, 1}, -472.5}});
  check_polynomial(6, 5, {{{1, 4, 1}, 51975.0}, {{3, 2, 1}, -103950.0}, {{5, 0, 1}, 10395.0}});
  check_polynomial(6, 6, {{{0, 6, 0}, -10395.0}, {{2, 4, 0}, 155925.0}, {{4, 2, 0}, -155925.0},
                          {{6, 0, 0}, 10395.0}});
}

TEST_CASE("real solid harmonics are harmonic", "[basis][solid_harmonics]") {
  for (int l = 0; l <= max_angular_momentum; ++l)
    for (int m = -l; m <= l; ++m) {
      INFO("(l, m) = (" << l << ", " << m << ")");
      const auto coefficients = real_solid_harmonic(l, m);
      double scale = 0.0;
      for (double c : coefficients) scale = std::max(scale, std::abs(c));
      REQUIRE(scale > 0.0);
      for (const auto& [powers, value] : laplacian(l, coefficients))
        CHECK_THAT(value, WithinAbs(0.0, 1e-12 * scale * l * l));
    }
  CHECK_THROWS_AS(real_solid_harmonic(2, 3), std::out_of_range);
  CHECK_THROWS_AS(real_solid_harmonic(2, -3), std::out_of_range);
}

TEST_CASE("order of the functions in a spherical shell", "[basis][solid_harmonics]") {
  CHECK(spherical_m_order(0) == std::vector<int>{0});
  CHECK(spherical_m_order(1) == std::vector<int>{1, -1, 0});  // x, y, z
  CHECK(spherical_m_order(2) == std::vector<int>{-2, -1, 0, 1, 2});
  CHECK(spherical_m_order(3) == std::vector<int>{-3, -2, -1, 0, 1, 2, 3});
}

TEST_CASE("spherical transformation gives orthonormal functions", "[basis][solid_harmonics]") {
  for (int l = 0; l <= max_angular_momentum; ++l) {
    INFO("l = " << l);
    const Eigen::MatrixXd& T = shell_transformation(l, true);
    const auto& powers = cartesian_powers(l);
    REQUIRE(T.rows() == n_spherical(l));
    REQUIRE(T.cols() == n_cartesian(l));

    // Overlap of the radially normalised cartesian functions on one centre:
    // integral of the product of the monomials over the unit sphere.
    Eigen::MatrixXd G(T.cols(), T.cols());
    for (Eigen::Index c = 0; c < T.cols(); ++c)
      for (Eigen::Index d = 0; d < T.cols(); ++d) {
        const auto& p = powers[static_cast<std::size_t>(c)];
        const auto& q = powers[static_cast<std::size_t>(d)];
        G(c, d) = sphere_monomial_integral(p.x + q.x, p.y + q.y, p.z + q.z);
      }
    const Eigen::MatrixXd overlap = T * G * T.transpose();
    CHECK((overlap - Eigen::MatrixXd::Identity(T.rows(), T.rows())).cwiseAbs().maxCoeff() < 1e-13);

    // Rows are positive multiples of the unnormalised harmonics.
    const auto m_order = spherical_m_order(l);
    for (Eigen::Index row = 0; row < T.rows(); ++row) {
      const auto polynomial = real_solid_harmonic(l, m_order[static_cast<std::size_t>(row)]);
      double ratio = 0.0;
      for (Eigen::Index c = 0; c < T.cols(); ++c) {
        const double a = polynomial[static_cast<std::size_t>(c)];
        if (std::abs(a) < 1e-12) {
          CHECK(T(row, c) == 0.0);
        } else if (ratio == 0.0) {
          ratio = T(row, c) / a;
          CHECK(ratio > 0.0);
        } else {
          CHECK_THAT(T(row, c) / a, WithinRel(ratio, 1e-12));
        }
      }
    }
  }
}

TEST_CASE("s and p transformations and the cartesian convention", "[basis][solid_harmonics]") {
  // s: 1/sqrt(4 pi); p: sqrt(3/(4 pi)) times x, y, z in this order.
  const Eigen::MatrixXd& s = shell_transformation(0, true);
  CHECK_THAT(s(0, 0), WithinRel(1.0 / std::sqrt(4.0 * M_PI), 1e-14));
  const Eigen::MatrixXd& p = shell_transformation(1, true);
  const Eigen::MatrixXd expected_p = std::sqrt(3.0 / (4.0 * M_PI)) * Eigen::MatrixXd::Identity(3, 3);
  CHECK((p - expected_p).cwiseAbs().maxCoeff() < 1e-14);

  // Cartesian shells: s and p as in the spherical case, identity from d on.
  CHECK(shell_transformation(0, false) == s);
  CHECK(shell_transformation(1, false) == p);
  for (int l = 2; l <= max_angular_momentum; ++l)
    CHECK(shell_transformation(l, false) == Eigen::MatrixXd::Identity(n_cartesian(l), n_cartesian(l)));
}
