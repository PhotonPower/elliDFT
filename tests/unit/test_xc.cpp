// Functionals: libxc behind the own interface, and the own LDA against it.
// Acceptance criterion M2: own LDA against libxc, energy density and
// derivatives agree to 1e-12.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "xc/functional.hpp"
#include "xc/xc_evaluator.hpp"

using namespace ellidft;
using Catch::Matchers::WithinRel;

namespace {

// Densities from 1e-10 to 1e4 electrons per cubic Bohr, logarithmically spaced.
std::vector<double> density_range() {
  std::vector<double> rho;
  for (int i = 0; i <= 280; ++i) rho.push_back(std::pow(10.0, -10.0 + 0.05 * i));
  return rho;
}

}  // namespace

TEST_CASE("functional data objects", "[xc]") {
  const Functional hf = Functional::by_name("hf");
  CHECK(hf.name == "HF");
  CHECK_FALSE(hf.has_density_part());
  CHECK(hf.exact_exchange == 1.0);

  const Functional lda = Functional::by_name("LDA");
  REQUIRE(lda.components.size() == 2);
  CHECK(lda.components[0].libxc_name == "LDA_X");
  CHECK(lda.components[1].libxc_name == "LDA_C_VWN");  // VWN5
  CHECK_FALSE(lda.needs_gradient);
  CHECK(lda.exact_exchange == 0.0);

  const Functional pbe = Functional::by_name("pbe");
  REQUIRE(pbe.components.size() == 2);
  CHECK(pbe.components[0].libxc_name == "GGA_X_PBE");
  CHECK(pbe.components[1].libxc_name == "GGA_C_PBE");
  CHECK(pbe.needs_gradient);
  CHECK(pbe.long_range_exchange == 0.0);

  CHECK_THROWS_AS(Functional::by_name("no-such-functional"), std::invalid_argument);
}

TEST_CASE("own LDA equals libxc: energy density and potential", "[xc][m2]") {
  const auto own = make_own_lda_evaluator();
  const auto reference = make_libxc_evaluator(Functional::by_name("LDA"));
  CHECK_FALSE(own->needs_gradient());
  CHECK_FALSE(reference->needs_gradient());

  const std::vector<double> rho = density_range();
  const std::size_t n = rho.size();
  std::vector<double> e1(n), v1(n), e2(n), v2(n);
  own->evaluate(n, rho.data(), nullptr, e1.data(), v1.data(), nullptr);
  reference->evaluate(n, rho.data(), nullptr, e2.data(), v2.data(), nullptr);

  double worst_energy = 0.0, worst_potential = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    // "Energiedichte": energy per particle, exc / rho.
    worst_energy = std::max(worst_energy, std::abs(e1[i] - e2[i]) / rho[i]);
    worst_potential = std::max(worst_potential, std::abs(v1[i] - v2[i]));
  }
  INFO("largest deviation: energy per particle " << worst_energy << ", potential " << worst_potential);
  CHECK(worst_energy <= 1e-12);
  CHECK(worst_potential <= 1e-12);
  SUCCEED("largest deviation: energy per particle " << worst_energy << ", potential " << worst_potential);

  // Exchange only, against the closed form -3/4 (3/pi)^(1/3) rho^(4/3).
  Functional exchange;
  exchange.components = {{"LDA_X", 1.0}};
  const auto slater = make_libxc_evaluator(exchange);
  slater->evaluate(n, rho.data(), nullptr, e2.data(), v2.data(), nullptr);
  for (std::size_t i = 0; i < n; ++i)
    CHECK_THAT(e2[i], WithinRel(-0.75 * std::cbrt(3.0 / M_PI) * std::pow(rho[i], 4.0 / 3.0), 1e-13));
}

TEST_CASE("potentials are derivatives of the energy density", "[xc]") {
  // d exc / d rho and d exc / d sigma by central differences. This fixes the
  // meaning of the outputs independently of the documentation of libxc.
  for (const char* name : {"LDA", "PBE"}) {
    INFO(name);
    const auto xc = make_libxc_evaluator(Functional::by_name(name));
    for (double rho : {1e-3, 0.05, 0.7, 12.0}) {
      for (double s : {0.0, 0.3, 2.0}) {  // sigma in units of rho^(8/3)
        const double sigma = s * std::pow(rho, 8.0 / 3.0);
        double e, v, vs, e_plus, e_minus, dummy1, dummy2;
        xc->evaluate(1, &rho, &sigma, &e, &v, &vs);

        const double h = 1e-4 * rho;
        const double rho_plus = rho + h, rho_minus = rho - h;
        xc->evaluate(1, &rho_plus, &sigma, &e_plus, &dummy1, &dummy2);
        xc->evaluate(1, &rho_minus, &sigma, &e_minus, &dummy1, &dummy2);
        INFO("rho = " << rho << ", sigma = " << sigma);
        CHECK_THAT((e_plus - e_minus) / (2.0 * h), WithinRel(v, 1e-6));

        if (xc->needs_gradient() && sigma > 0.0) {
          const double hs = 1e-4 * sigma;
          const double sigma_plus = sigma + hs, sigma_minus = sigma - hs;
          xc->evaluate(1, &rho, &sigma_plus, &e_plus, &dummy1, &dummy2);
          xc->evaluate(1, &rho, &sigma_minus, &e_minus, &dummy1, &dummy2);
          // Truncation error of the central difference (relative 1e-6) plus
          // its rounding error: the sigma-dependent part of the energy can be
          // a millionth of the energy itself, and e is known to about 2 ulp.
          const double tolerance = 1e-6 * std::abs(vs) + 4e-16 * std::abs(e) / hs;
          CHECK_THAT((e_plus - e_minus) / (2.0 * hs), Catch::Matchers::WithinAbs(vs, tolerance));
        } else {
          CHECK((xc->needs_gradient() || vs == 0.0));
        }
      }
    }
  }
  // The own LDA likewise.
  const auto own = make_own_lda_evaluator();
  for (double rho : {1e-3, 0.05, 0.7, 12.0}) {
    double e, v, e_plus, e_minus, dummy;
    own->evaluate(1, &rho, nullptr, &e, &v, nullptr);
    const double h = 1e-4 * rho, rho_plus = rho + h, rho_minus = rho - h;
    own->evaluate(1, &rho_plus, nullptr, &e_plus, &dummy, nullptr);
    own->evaluate(1, &rho_minus, nullptr, &e_minus, &dummy, nullptr);
    CHECK_THAT((e_plus - e_minus) / (2.0 * h), WithinRel(v, 1e-6));
  }
}

TEST_CASE("functional evaluators: limits and errors", "[xc]") {
  const auto own = make_own_lda_evaluator();
  const double zero = 0.0;
  double e = 1.0, v = 1.0, vs = 1.0;
  own->evaluate(1, &zero, nullptr, &e, &v, &vs);
  CHECK(e == 0.0);
  CHECK(v == 0.0);
  CHECK(vs == 0.0);

  // A coefficient scales the result.
  Functional half = Functional::by_name("LDA");
  for (auto& c : half.components) c.coefficient = 0.5;
  const auto full = make_libxc_evaluator(Functional::by_name("LDA"));
  const auto scaled = make_libxc_evaluator(half);
  const double rho = 0.3;
  double e_full, v_full, e_half, v_half;
  full->evaluate(1, &rho, nullptr, &e_full, &v_full, nullptr);
  scaled->evaluate(1, &rho, nullptr, &e_half, &v_half, nullptr);
  CHECK_THAT(e_half, WithinRel(0.5 * e_full, 1e-14));
  CHECK_THAT(v_half, WithinRel(0.5 * v_full, 1e-14));

  Functional unknown;
  unknown.components = {{"NOT_A_FUNCTIONAL", 1.0}};
  CHECK_THROWS_AS(make_libxc_evaluator(unknown), std::invalid_argument);
  Functional meta;
  meta.components = {{"MGGA_X_TPSS", 1.0}};
  CHECK_THROWS_AS(make_libxc_evaluator(meta), std::invalid_argument);
  // GGA without sigma
  const auto pbe = make_libxc_evaluator(Functional::by_name("PBE"));
  CHECK_THROWS_AS(pbe->evaluate(1, &rho, nullptr, &e, &v, nullptr), std::invalid_argument);
  CHECK(std::string(libxc_version()).size() > 0);
}
