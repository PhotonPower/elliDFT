// Absorbing oscillator models and effective media.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>

#include "ellipsometry/materials.hpp"

using namespace ellidft::ellipsometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Lorentz oscillators", "[ellipsometry]") {
  LorentzMaterial m(2.0, {{12.0, 4.0, 0.5}, {3.0, 6.5, 0.2}});
  // eps = eps_inf + sum s / (E0^2 - E^2 - i Gamma E)
  const double E = 3.0;
  const Complex expected = 2.0 + 12.0 / Complex(16.0 - 9.0, -1.5) + 3.0 / Complex(42.25 - 9.0, -0.6);
  CHECK_THAT(std::abs(m.permittivity(E) - expected), WithinAbs(0.0, 1e-14));
  CHECK(m.permittivity(E).imag() > 0.0);                     // absorption, never gain
  CHECK(m.index(1239.841984 / E).imag() > 0.0);
  CHECK_THAT(std::abs(m.index(1239.841984 / E) * m.index(1239.841984 / E) - expected), WithinAbs(0.0, 1e-13));
  // static limit: eps_inf + sum s / E0^2
  CHECK_THAT(m.permittivity(0.0).real(), WithinRel(2.0 + 12.0 / 16.0 + 3.0 / 42.25, 1e-14));
  CHECK(photon_energy_ev(620.0) == 1239.841984 / 620.0);
  // parameters
  CHECK(m.parameter_names() == std::vector<std::string>{"eps_inf", "strength1", "energy1", "width1", "strength2", "energy2", "width2"});
  CHECK(m.parameter(5) == 6.5);
  m.set_parameter(3, 0.7);
  CHECK(m.parameter(3) == 0.7);
}

TEST_CASE("Kramers-Kronig integration reproduces the real part of a Lorentz oscillator", "[ellipsometry]") {
  // eps2 of one oscillator is known in closed form together with its real part.
  const double s = 9.0, E0 = 4.2, G = 0.6;
  const LorentzMaterial m(1.0, {{s, E0, G}});
  auto eps2 = [&](double E) { return m.permittivity(E).imag(); };
  double largest = 0.0;
  for (double E : {0.0, 1.0, 2.5, 4.0, 4.2, 4.5, 6.0, 10.0}) {
    const double numerical = kramers_kronig_real_part(eps2, 0.0, E0, E);
    largest = std::max(largest, std::abs(numerical - (m.permittivity(E).real() - 1.0)));
  }
  std::printf("Kramers-Kronig integration, Lorentz oscillator: largest deviation of the real part %.1e\n", largest);
  CHECK(largest < 1e-6);
}

namespace {

// Closed form of the real part of the Tauc-Lorentz model as given on the Wikipedia page (after Jellison and
// Modine 1996). The page's formulas came through the retrieval tool with lost radicals and exponents and were
// restored by it; hence only a check here, not part of the library. In the last logarithm the retrieved text
// had no square root over (E0^2 - Eg^2)^2 + Eg^2 C^2; with the root the formula agrees with the numerical
// integration to 1e-11, without it it does not.
double tauc_lorentz_eps1_closed(double eps_inf, double A, double E0, double C, double Eg, double E) {
  const double pi = M_PI;
  const double a_ln = (Eg * Eg - E0 * E0) * E * E + Eg * Eg * C * C - E0 * E0 * (E0 * E0 + 3.0 * Eg * Eg);
  const double a_atan = (E * E - E0 * E0) * (E0 * E0 + Eg * Eg) + Eg * Eg * C * C;
  const double alpha = std::sqrt(4.0 * E0 * E0 - C * C), gamma = std::sqrt(E0 * E0 - C * C / 2.0);
  const double zeta4 = std::pow(E * E - gamma * gamma, 2) + alpha * alpha * C * C / 4.0;
  return eps_inf +
         A * C * a_ln / (2.0 * pi * zeta4 * alpha * E0) * std::log((E0 * E0 + Eg * Eg + alpha * Eg) / (E0 * E0 + Eg * Eg - alpha * Eg)) -
         A * a_atan / (pi * zeta4 * E0) * (pi - std::atan((alpha + 2.0 * Eg) / C) + std::atan((alpha - 2.0 * Eg) / C)) +
         2.0 * A * E0 * Eg / (pi * zeta4 * alpha) * (E * E - gamma * gamma) * (pi + 2.0 * std::atan(2.0 * (gamma * gamma - Eg * Eg) / (alpha * C))) -
         A * E0 * C / (pi * zeta4) * (E * E + Eg * Eg) / E * std::log(std::abs(E - Eg) / (E + Eg)) +
         2.0 * A * E0 * C * Eg / (pi * zeta4) * std::log(std::abs(E - Eg) * (E + Eg) / std::sqrt(std::pow(E0 * E0 - Eg * Eg, 2) + Eg * Eg * C * C));
}

}  // namespace

TEST_CASE("Tauc-Lorentz model", "[ellipsometry]") {
  const double eps_inf = 1.5, A = 120.0, E0 = 3.8, C = 1.4, Eg = 1.6;
  TaucLorentzMaterial m(eps_inf, A, E0, C, Eg);
  // no absorption below the gap, continuous at the gap, positive above
  CHECK(m.eps2(1.0) == 0.0);
  CHECK(m.eps2(Eg) == 0.0);
  CHECK(m.eps2(Eg + 1e-6) < 1e-9);
  CHECK(m.eps2(3.0) > 0.0);
  const double E = 3.0, d = E * E - E0 * E0;
  CHECK_THAT(m.eps2(E), WithinRel(A * E0 * C * (E - Eg) * (E - Eg) / (E * (d * d + C * C * E * E)), 1e-14));
  CHECK(m.index(1239.841984 / 1.0).imag() == 0.0);   // transparent below the gap
  CHECK(m.index(1239.841984 / 3.0).imag() > 0.0);
  // the real part exceeds eps_inf below the absorption and tends to eps_inf far above it
  CHECK(m.permittivity(0.5).real() > eps_inf);
  CHECK_THAT(m.permittivity(200.0).real(), WithinAbs(eps_inf, 0.02));
  // numerical Kramers-Kronig integration against the closed form of the source
  double largest = 0.0;
  for (double energy : {0.5, 1.0, 1.5, 2.0, 3.0, 3.8, 5.0, 8.0}) {
    const double numerical = m.permittivity(energy).real();
    const double closed = tauc_lorentz_eps1_closed(eps_inf, A, E0, C, Eg, energy);
    largest = std::max(largest, std::abs(numerical - closed));
    INFO("E = " << energy << ": numerical " << numerical << ", closed form " << closed);
    CHECK_THAT(numerical, WithinAbs(closed, 2e-4));
  }
  std::printf("Tauc-Lorentz: largest deviation between numerical integration and the closed form %.1e\n", largest);
  CHECK(m.parameter_names().size() == 5);
  m.set_parameter(4, 2.0);
  CHECK(m.eps2(1.9) == 0.0);
  m.set_parameter(3, 0.0);
  CHECK_THROWS_AS(m.index(500.0), std::domain_error);
}

TEST_CASE("Bruggeman effective medium", "[ellipsometry]") {
  const Complex a(2.25, 0.0), b(1.0, 0.0), absorbing(-8.0, 3.0);
  // limits
  CHECK_THAT(std::abs(BruggemanMaterial::mix(a, b, 0.0) - a), WithinAbs(0.0, 1e-14));
  CHECK_THAT(std::abs(BruggemanMaterial::mix(a, b, 1.0) - b), WithinAbs(0.0, 1e-14));
  CHECK_THAT(std::abs(BruggemanMaterial::mix(a, a, 0.37) - a), WithinAbs(0.0, 1e-14));
  for (const Complex other : {b, absorbing, Complex(4.0, 0.5)}) {
    for (double f : {0.1, 0.33, 0.5, 0.8}) {
      const Complex eps = BruggemanMaterial::mix(a, other, f);
      // defining equation: sum of f_i (eps_i - eps) / (eps_i + 2 eps) = 0
      const Complex residual = (1.0 - f) * (a - eps) / (a + 2.0 * eps) + f * (other - eps) / (other + 2.0 * eps);
      CHECK_THAT(std::abs(residual), WithinAbs(0.0, 1e-13));
      CHECK(eps.imag() >= -1e-13);                                   // never gain
      // the same mixture with the components exchanged
      CHECK_THAT(std::abs(BruggemanMaterial::mix(other, a, 1.0 - f) - eps), WithinAbs(0.0, 1e-13));
    }
  }
  // transparent components: between the two
  const Complex half = BruggemanMaterial::mix(a, b, 0.5);
  CHECK(half.real() > 1.0);
  CHECK(half.real() < 2.25);
  CHECK(half.imag() == 0.0);
  // dilute limit: d eps / d f = 3 eps_a (eps_b - eps_a) / (eps_b + 2 eps_a)
  const double f = 1e-6;
  CHECK_THAT((BruggemanMaterial::mix(a, b, f) - a).real() / f, WithinRel((3.0 * a * (b - a) / (b + 2.0 * a)).real(), 1e-4));

  // as a material: roughness as half film, half ambient
  auto film = std::make_shared<ConstantMaterial>(1.5, 0.0);
  auto air = std::make_shared<ConstantMaterial>(1.0, 0.0);
  BruggemanMaterial rough(film, air, 0.5);
  CHECK_THAT(std::norm(rough.index(500.0)), WithinRel(half.real(), 1e-14));
  CHECK(rough.parameter_names() == std::vector<std::string>{"fraction_b"});
  rough.set_parameter(0, 0.0);
  CHECK_THAT(rough.index(500.0).real(), WithinRel(1.5, 1e-14));
  film->set_parameter(0, 1.6);   // the components stay adjustable through their own objects
  CHECK_THAT(rough.index(500.0).real(), WithinRel(1.6, 1e-14));
  CHECK_THROWS_AS(BruggemanMaterial::mix(a, b, 1.5), std::domain_error);
  CHECK_THROWS_AS(BruggemanMaterial(nullptr, air, 0.5), std::invalid_argument);
}
