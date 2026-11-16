// Dispersion models fitted to a refractive index: known coefficients must come back.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "ellipsometry/dispersion_fit.hpp"
#include "ellipsometry/materials.hpp"

using namespace ellidft::ellipsometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {
const std::vector<double> wavelengths = {400, 450, 500, 550, 589, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
}

TEST_CASE("dispersion fit: Sellmeier coefficients of polystyrene and PMMA are recovered", "[dispersion_fit]") {
  // One-term Sellmeier fits of the refractiveindex.info database (docs/sources/README.md), as in
  // tests/reference/gen_m6.py: polystyrene B = 1.4435, C = 0.020216; PMMA B = 1.1819, C = 0.011313.
  for (const auto& [B, C] : {std::pair{1.4435, 0.020216}, std::pair{1.1819, 0.011313}}) {
    const SellmeierMaterial material({B}, {C});
    std::vector<double> n;
    for (double w : wavelengths) n.push_back(material.index(w).real());
    const SellmeierFit fit = fit_sellmeier(wavelengths, n);
    CHECK(fit.converged);
    CHECK_THAT(fit.B, WithinRel(B, 1e-7));
    CHECK_THAT(fit.C_um2, WithinRel(C, 1e-6));
    CHECK(fit.largest < 1e-9);
    // the Cauchy form describes the same data to a few 1e-5
    const CauchyFit cauchy = fit_cauchy(wavelengths, n);
    CHECK(cauchy.largest < 1e-4);
    CHECK(cauchy.A > 1.4);
    CHECK(cauchy.B > 0.0);
  }
}

TEST_CASE("dispersion fit: Cauchy coefficients are recovered", "[dispersion_fit]") {
  const CauchyMaterial material(1.5632, 0.0092, 0.0002);
  std::vector<double> n;
  for (double w : wavelengths) n.push_back(material.index(w).real());
  const CauchyFit fit = fit_cauchy(wavelengths, n);
  CHECK_THAT(fit.A, WithinAbs(1.5632, 1e-10));
  CHECK_THAT(fit.B, WithinAbs(0.0092, 1e-10));
  CHECK_THAT(fit.C, WithinAbs(0.0002, 1e-10));
  CHECK(fit.rms < 1e-12);
  // two or three points: no fourth-power term
  const CauchyFit two = fit_cauchy({500.0, 700.0}, {1.60, 1.58});
  CHECK(two.C == 0.0);
  CHECK_THAT(two.A + two.B / 0.25, WithinAbs(1.60, 1e-12));
  CHECK_THAT(two.A + two.B / 0.49, WithinAbs(1.58, 1e-12));
}

TEST_CASE("dispersion fit: invalid input", "[dispersion_fit]") {
  CHECK_THROWS_AS(fit_cauchy({500.0}, {1.5}), std::invalid_argument);
  CHECK_THROWS_AS(fit_cauchy({500.0, 600.0}, {1.5}), std::invalid_argument);
  CHECK_THROWS_AS(fit_cauchy({500.0, -600.0}, {1.5, 1.5}), std::invalid_argument);
  CHECK_THROWS_AS(fit_sellmeier({500.0, 600.0}, {1.5, 1.5}), std::invalid_argument);
  CHECK_THROWS_AS(fit_sellmeier({500.0, 600.0, 700.0}, {1.5, 1.0, 1.4}), std::invalid_argument);
}
