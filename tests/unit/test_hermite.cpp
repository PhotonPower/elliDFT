#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "integrals/boys.hpp"
#include "integrals/hermite.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Hermite Gaussians Lambda_0..Lambda_tmax at x - P, from
//   x_P Lambda_t = Lambda_{t+1} / (2p) + t Lambda_{t-1}
// (Hernández Vera 2025, Eq. 14) and Lambda_0 = exp(-p x_P^2) (Eq. 13).
std::vector<double> hermite_gaussians(int tmax, double p, double xP) {
  std::vector<double> lambda(static_cast<std::size_t>(tmax + 1));
  lambda[0] = std::exp(-p * xP * xP);
  for (int t = 0; t < tmax; ++t) {
    const double previous = t > 0 ? lambda[static_cast<std::size_t>(t - 1)] : 0.0;
    lambda[static_cast<std::size_t>(t + 1)] = 2.0 * p * (xP * lambda[static_cast<std::size_t>(t)] - t * previous);
  }
  return lambda;
}

template <typename Function>
double simpson(Function&& f, double lower, double upper, int intervals) {
  const double h = (upper - lower) / intervals;
  double sum = f(lower) + f(upper);
  for (int i = 1; i < intervals; ++i) sum += (i % 2 == 1 ? 4.0 : 2.0) * f(lower + i * h);
  return sum * h / 3.0;
}

}  // namespace

TEST_CASE("Hermite coefficients reproduce the product of two Gaussians", "[integrals][hermite]") {
  const int imax = 6, jmax = 8;
  struct Case { double a, b, A, B; };
  for (const Case& c : {Case{0.9, 1.7, -0.3, 0.8}, Case{2.5, 0.2, 1.1, -0.4}, Case{0.6, 0.6, 0.5, 0.5}}) {
    const HermiteE E(imax, jmax, c.a, c.b, c.A, c.B);
    const double p = c.a + c.b;
    const double P = (c.a * c.A + c.b * c.B) / p;
    for (double x : {-1.7, -0.2, 0.33, 0.9, 2.1}) {
      const auto lambda = hermite_gaussians(imax + jmax, p, x - P);
      for (int i = 0; i <= imax; ++i)
        for (int j = 0; j <= jmax; ++j) {
          const double product = std::pow(x - c.A, i) * std::pow(x - c.B, j) *
                                 std::exp(-c.a * (x - c.A) * (x - c.A) - c.b * (x - c.B) * (x - c.B));
          double expansion = 0.0, magnitude = 0.0;
          for (int t = 0; t <= i + j; ++t) {
            expansion += E(i, j, t) * lambda[static_cast<std::size_t>(t)];
            magnitude += std::abs(E(i, j, t) * lambda[static_cast<std::size_t>(t)]);
          }
          INFO("a=" << c.a << " b=" << c.b << " x=" << x << " i=" << i << " j=" << j);
          CHECK_THAT(expansion, WithinAbs(product, 1e-12 * magnitude + 1e-300));
        }
    }
  }
}

TEST_CASE("Hermite coefficients: boundary values and symmetry", "[integrals][hermite]") {
  const double a = 1.3, b = 0.7, A = 0.2, B = -0.9;
  const HermiteE E(4, 4, a, b, A, B);
  // Eq. 19
  CHECK_THAT(E(0, 0, 0), WithinRel(std::exp(-a * b / (a + b) * (A - B) * (A - B)), 1e-15));
  // Eq. 20
  CHECK(E(2, 3, -1) == 0.0);
  CHECK(E(2, 3, 6) == 0.0);
  // Highest coefficient: E_{i+j}^{ij} = E_0^{00} / (2p)^(i+j), from Eqs. 17, 18.
  for (int i = 0; i <= 4; ++i)
    for (int j = 0; j <= 4; ++j)
      CHECK_THAT(E(i, j, i + j), WithinRel(E(0, 0, 0) / std::pow(2.0 * (a + b), i + j), 1e-13));
  // Exchanging the two functions exchanges the indices. A coefficient is a
  // sum of terms of both signs (P-A and P-B have opposite signs), so the
  // rounding error scales with the largest coefficient, not with the value.
  const HermiteE swapped(4, 4, b, a, B, A);
  double scale = 0.0;
  for (int i = 0; i <= 4; ++i)
    for (int j = 0; j <= 4; ++j)
      for (int t = 0; t <= i + j; ++t) scale = std::max(scale, std::abs(E(i, j, t)));
  for (int i = 0; i <= 4; ++i)
    for (int j = 0; j <= 4; ++j)
      for (int t = 0; t <= i + j; ++t)
        CHECK_THAT(swapped(j, i, t), WithinAbs(E(i, j, t), 1e-14 * scale));
}

// Independent check of the one-dimensional overlap (Hernández Vera 2025,
// Eq. 35) by numerical quadrature.
TEST_CASE("one-dimensional overlap E_0 sqrt(pi/p) against quadrature", "[integrals][hermite]") {
  const double a = 0.8, b = 1.4, A = -0.6, B = 0.7;
  const int imax = 6, jmax = 8;
  const HermiteE E(imax, jmax, a, b, A, B);
  const double p = a + b;
  for (int i = 0; i <= imax; ++i)
    for (int j = 0; j <= jmax; ++j) {
      const double numeric = simpson(
          [&](double x) {
            return std::pow(x - A, i) * std::pow(x - B, j) *
                   std::exp(-a * (x - A) * (x - A) - b * (x - B) * (x - B));
          },
          -12.0, 12.0, 40000);
      INFO("i=" << i << " j=" << j);
      // Integrals of odd integrands are small differences of large parts; the
      // absolute term covers them.
      CHECK_THAT(E(i, j, 0) * std::sqrt(M_PI / p), WithinAbs(numeric, 1e-10 * std::abs(numeric) + 1e-12));
    }
}

TEST_CASE("Hermite index list", "[integrals][hermite]") {
  for (int L = 0; L <= 12; ++L) {
    const auto list = hermite_indices(L);
    REQUIRE(list.size() == static_cast<std::size_t>(n_hermite(L)));
    for (const auto& h : list) {
      CHECK(h.t >= 0);
      CHECK(h.u >= 0);
      CHECK(h.v >= 0);
      CHECK(h.t + h.u + h.v <= L);
    }
  }
  CHECK(n_hermite(0) == 1);
  CHECK(n_hermite(1) == 4);
  CHECK(n_hermite(2) == 10);
}

TEST_CASE("Hermite Coulomb integrals: source term and derivatives", "[integrals][hermite]") {
  const int L = 8;
  const double alpha = 0.83;
  const double X = 0.41, Y = -0.77, Z = 0.29;
  HermiteCoulomb R(L), plus(L), minus(L);
  R.compute(L, alpha, X, Y, Z);

  SECTION("R_000 is the Boys function F_0(alpha R^2)") {
    double f0;
    boys(0, alpha * (X * X + Y * Y + Z * Z), &f0);
    CHECK_THAT(R.value(0, 0, 0), WithinRel(f0, 1e-15));
  }
  SECTION("a lower L gives the same values") {
    HermiteCoulomb small(L);
    small.compute(3, alpha, X, Y, Z);
    for (const auto& h : hermite_indices(3))
      CHECK_THAT(small.value(h.t, h.u, h.v), WithinRel(R.value(h.t, h.u, h.v), 1e-13));
  }
  // Definition (Hernández Vera 2025, Eq. 50): R_{t+1,u,v} is the derivative of
  // R_{tuv} with respect to the x component of the distance vector, and the
  // same for y and z. Central differences with Richardson extrapolation.
  SECTION("R_{tuv} are successive derivatives") {
    auto derivative = [&](int axis, int t, int u, int v) {
      auto difference = [&](double h) {
        const double dx = axis == 0 ? h : 0.0, dy = axis == 1 ? h : 0.0, dz = axis == 2 ? h : 0.0;
        plus.compute(L - 1, alpha, X + dx, Y + dy, Z + dz);
        minus.compute(L - 1, alpha, X - dx, Y - dy, Z - dz);
        return (plus.value(t, u, v) - minus.value(t, u, v)) / (2.0 * h);
      };
      const double h = 2e-3;
      return (4.0 * difference(h / 2.0) - difference(h)) / 3.0;
    };
    double scale = 0.0;
    for (const auto& h : hermite_indices(L)) scale = std::max(scale, std::abs(R.value(h.t, h.u, h.v)));
    for (const auto& h : hermite_indices(L - 1)) {
      INFO("t=" << h.t << " u=" << h.u << " v=" << h.v);
      const double tolerance = 1e-8 * scale;
      CHECK_THAT(derivative(0, h.t, h.u, h.v), WithinAbs(R.value(h.t + 1, h.u, h.v), tolerance));
      CHECK_THAT(derivative(1, h.t, h.u, h.v), WithinAbs(R.value(h.t, h.u + 1, h.v), tolerance));
      CHECK_THAT(derivative(2, h.t, h.u, h.v), WithinAbs(R.value(h.t, h.u, h.v + 1), tolerance));
    }
  }
  SECTION("zero distance: odd orders vanish") {
    R.compute(L, alpha, 0.0, 0.0, 0.0);
    CHECK(R.value(0, 0, 0) == 1.0);
    for (const auto& h : hermite_indices(L))
      if (h.t % 2 == 1 || h.u % 2 == 1 || h.v % 2 == 1) CHECK(R.value(h.t, h.u, h.v) == 0.0);
    // Second derivative of F_0(alpha x^2) at 0: -2 alpha F_1(0) = -2 alpha / 3.
    CHECK_THAT(R.value(2, 0, 0), WithinRel(-2.0 * alpha / 3.0, 1e-14));
    CHECK_THAT(R.value(0, 0, 2), WithinRel(-2.0 * alpha / 3.0, 1e-14));
  }
  SECTION("unsupported total angular momentum") {
    CHECK_THROWS_AS(HermiteCoulomb(boys_max_order + 1), std::out_of_range);
    CHECK_THROWS_AS(HermiteCoulomb(-1), std::out_of_range);
  }
}
