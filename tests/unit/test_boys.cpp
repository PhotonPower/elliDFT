// Boys function against the high-precision reference
// (tests/reference/data/m1/boys.json, 60-digit decimal arithmetic).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>

#include <nlohmann/json.hpp>

#include "integrals/boys.hpp"

using namespace ellidft;

namespace {

struct BoysReference {
  std::vector<double> x;
  std::vector<std::vector<double>> values;  // values[i][n] = F_n(x[i])
};

const BoysReference& reference() {
  static const BoysReference ref = [] {
    std::ifstream in(std::string(ELLIDFT_SOURCE_DIR) + "/tests/reference/data/m1/boys.json");
    REQUIRE(in.good());
    const nlohmann::json j = nlohmann::json::parse(in);
    return BoysReference{j.at("x").get<std::vector<double>>(),
                         j.at("values").get<std::vector<std::vector<double>>>()};
  }();
  return ref;
}

// Largest relative error over all reference points with x in [x_min, x_max]
// and all orders up to n_max. `call_order` is the nmax passed to the function
// under test (it influences the result through the recursions).
template <typename Function>
double max_relative_error(Function&& f, int call_order, int n_max, double x_min, double x_max) {
  const BoysReference& ref = reference();
  double worst = 0.0;
  std::vector<double> out(static_cast<std::size_t>(call_order + 1));
  for (std::size_t i = 0; i < ref.x.size(); ++i) {
    if (ref.x[i] < x_min || ref.x[i] > x_max) continue;
    f(call_order, ref.x[i], out.data());
    for (int n = 0; n <= n_max; ++n) {
      const double exact = ref.values[i][static_cast<std::size_t>(n)];
      worst = std::max(worst, std::abs(out[static_cast<std::size_t>(n)] - exact) / exact);
    }
  }
  return worst;
}

}  // namespace

// Acceptance criterion M1: F_n(x) for n <= 24 and 0 <= x <= 500, relative
// error at most 1e-13 against a calculation with high precision.
TEST_CASE("Boys function: acceptance range n <= 24, 0 <= x <= 500", "[boys][m1]") {
  const double worst = max_relative_error(boys, 24, 24, 0.0, 500.0);
  INFO("largest relative error: " << worst);
  CHECK(worst <= 1e-13);
}

TEST_CASE("Boys function: all supported orders and very large arguments", "[boys]") {
  REQUIRE(reference().values.front().size() == static_cast<std::size_t>(boys_max_order + 1));
  const double worst = max_relative_error(boys, boys_max_order, boys_max_order, 0.0, 1.0e5);
  INFO("largest relative error: " << worst);
  CHECK(worst <= 1e-13);
}

// The result for order n must not depend on how many orders are requested
// beyond rounding: every nmax takes a different path through the recursions.
TEST_CASE("Boys function: every requested maximum order", "[boys]") {
  for (int nmax = 0; nmax <= boys_max_order; ++nmax) {
    const double worst = max_relative_error(boys, nmax, nmax, 0.0, 1.0e5);
    INFO("nmax = " << nmax << ", largest relative error: " << worst);
    CHECK(worst <= 1e-13);
  }
}

TEST_CASE("Boys function: series implementation used for the table", "[boys]") {
  // The table needs orders up to boys_max_order + 7; the reference stops at
  // boys_max_order, so the higher orders are checked through the downward
  // recursion they feed.
  const double worst = max_relative_error(boys_series, boys_max_order, boys_max_order, 0.0, 40.0);
  INFO("largest relative error: " << worst);
  CHECK(worst <= 2e-15);
}

TEST_CASE("Boys function: exact values and analytic properties", "[boys]") {
  double f[boys_max_order + 1];

  SECTION("F_n(0) = 1/(2n+1)") {
    boys(boys_max_order, 0.0, f);
    for (int n = 0; n <= boys_max_order; ++n) CHECK(f[n] == 1.0 / (2 * n + 1));
  }
  SECTION("F_0(x) = sqrt(pi/(4x)) erf(sqrt(x))") {
    for (double x : {1e-3, 0.05, 0.7, 3.3, 12.0, 35.95, 36.0, 36.05, 80.0, 499.0}) {
      boys(0, x, f);
      const double exact = std::sqrt(M_PI / (4.0 * x)) * std::erf(std::sqrt(x));
      INFO("x = " << x);
      CHECK(std::abs(f[0] - exact) <= 4e-15 * exact);
    }
  }
  SECTION("continuity across table points and the switch to the asymptotic branch") {
    // Neighbouring arguments on both sides of a method boundary: the function
    // values may differ only by the slope, |dF_n/dx| = F_(n+1) <= F_n.
    for (double x0 : {0.05, 0.15, 17.95, 35.95, 36.0}) {
      double below[boys_max_order + 1], above[boys_max_order + 1];
      const double h = 1e-9;
      boys(boys_max_order, x0 - h, below);
      boys(boys_max_order, x0 + h, above);
      for (int n = 0; n <= boys_max_order; ++n) {
        INFO("x0 = " << x0 << ", n = " << n);
        CHECK(std::abs(above[n] - below[n]) <= (2.0 * h + 2e-13) * below[n]);
      }
    }
  }
  SECTION("0 < F_(n+1)(x) < F_n(x) <= 1/(2n+1)") {
    for (double x = 0.0; x < 120.0; x += 0.731) {
      boys(boys_max_order, x, f);
      for (int n = 0; n <= boys_max_order; ++n) {
        CHECK(f[n] > 0.0);
        CHECK(f[n] <= 1.0 / (2 * n + 1));
        if (n > 0) CHECK(f[n] < f[n - 1]);
      }
    }
  }
}
