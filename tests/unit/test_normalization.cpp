#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <sstream>
#include <vector>

#include "basis/basis_set.hpp"

using namespace ellidft;
using Catch::Matchers::WithinRel;

namespace {

// Composite Simpson rule on [0, upper]; independent of any closed formula.
template <typename Function>
double simpson(Function&& f, double upper, int intervals) {
  const double h = upper / intervals;
  double sum = f(0.0) + f(upper);
  for (int i = 1; i < intervals; ++i) sum += (i % 2 == 1 ? 4.0 : 2.0) * f(i * h);
  return sum * h / 3.0;
}

double radial_norm(int l, const std::vector<double>& exponents, const std::vector<double>& d) {
  double smallest = exponents[0];
  for (double a : exponents) smallest = std::min(smallest, a);
  // exp(-2 a r^2) r^(2l+2) is below 1e-30 of its maximum well before this bound
  const double upper = 14.0 / std::sqrt(smallest);
  return simpson(
      [&](double r) {
        double R = 0.0;
        for (std::size_t k = 0; k < d.size(); ++k) R += d[k] * std::exp(-exponents[k] * r * r);
        return R * R * std::pow(r, 2 * l + 2);
      },
      upper, 200000);
}

}  // namespace

TEST_CASE("radial Gaussian integral against numerical quadrature", "[basis][normalization]") {
  for (int l = 0; l <= max_angular_momentum; ++l) {
    for (double zeta : {0.05, 0.8, 3.7, 150.0}) {
      INFO("l = " << l << ", zeta = " << zeta);
      const double numeric = simpson(
          [&](double r) { return std::pow(r, 2 * l + 2) * std::exp(-zeta * r * r); },
          12.0 / std::sqrt(zeta), 200000);
      CHECK_THAT(radial_gaussian_integral(l, zeta), WithinRel(numeric, 1e-11));
    }
  }
  // Closed values: l = 0 gives sqrt(pi)/(4 zeta^(3/2)).
  CHECK_THAT(radial_gaussian_integral(0, 1.0), WithinRel(std::sqrt(M_PI) / 4.0, 1e-15));
  CHECK_THAT(radial_gaussian_integral(0, 4.0), WithinRel(std::sqrt(M_PI) / 32.0, 1e-15));
}

TEST_CASE("normalised contractions have unit radial norm", "[basis][normalization]") {
  SECTION("single primitives, all angular momenta") {
    for (int l = 0; l <= max_angular_momentum; ++l) {
      for (double a : {0.04, 0.9, 25.0}) {
        INFO("l = " << l << ", exponent = " << a);
        // The coefficient from the file must not matter for a single primitive,
        // apart from its sign.
        const auto d = normalized_contraction(l, {a}, {0.37});
        const auto d_negative = normalized_contraction(l, {a}, {-2.0});
        CHECK_THAT(radial_norm(l, {a}, d), WithinRel(1.0, 1e-10));
        CHECK_THAT(d_negative[0], WithinRel(-d[0], 1e-14));
        CHECK(d[0] > 0.0);
      }
    }
  }
  SECTION("contractions") {
    const std::vector<double> exponents = {130.7, 23.8, 6.44, 0.52};
    const std::vector<double> coefficients = {0.019, 0.138, 0.478, -0.501};
    for (int l = 0; l <= max_angular_momentum; ++l) {
      INFO("l = " << l);
      const auto d = normalized_contraction(l, exponents, coefficients);
      CHECK_THAT(radial_norm(l, exponents, d), WithinRel(1.0, 1e-10));
      // Scaling all file coefficients changes nothing.
      std::vector<double> scaled = coefficients;
      for (double& c : scaled) c *= 7.3;
      const auto d_scaled = normalized_contraction(l, exponents, scaled);
      for (std::size_t k = 0; k < d.size(); ++k) CHECK_THAT(d_scaled[k], WithinRel(d[k], 1e-13));
    }
  }
  SECTION("s primitive: closed form (2a/pi)^(3/4) * sqrt(4 pi)") {
    // A normalised s Gaussian is (2a/pi)^(3/4) exp(-a r^2); the angular factor
    // 1/sqrt(4 pi) is not part of the radial function.
    const double a = 1.3;
    const auto d = normalized_contraction(0, {a}, {1.0});
    CHECK_THAT(d[0], WithinRel(std::pow(2.0 * a / M_PI, 0.75) * std::sqrt(4.0 * M_PI), 1e-14));
  }
  SECTION("invalid input") {
    CHECK_THROWS_AS(normalized_contraction(0, {1.0, 2.0}, {1.0}), std::invalid_argument);
    CHECK_THROWS_AS(normalized_contraction(0, {1.0}, {0.0}), std::invalid_argument);
  }
}

TEST_CASE("basis of a molecule carries normalised coefficients", "[basis][normalization]") {
  std::istringstream in(
      "BASIS \"ao basis\" SPHERICAL\nH S\n 3.0 0.2\n 0.5 0.8\nH D\n 1.1 1.0\nEND\n");
  const BasisLibrary library = BasisLibrary::parse_nwchem(in, "test");
  Molecule mol;
  mol.atoms = {{1, {0.0, 0.0, 0.0}}};
  const BasisSet basis = BasisSet::build(mol, library);
  REQUIRE(basis.n_shells() == 2);
  for (const Shell& shell : basis.shells()) {
    REQUIRE(shell.norm_coefficients.size() == shell.n_primitives());
    CHECK(shell.norm_coefficients == normalized_contraction(shell.l, shell.exponents, shell.coefficients));
    CHECK_THAT(radial_norm(shell.l, shell.exponents, shell.norm_coefficients), WithinRel(1.0, 1e-10));
  }
}
