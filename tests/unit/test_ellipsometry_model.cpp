// Module ellipsometry: material models and the forward model (reflection of
// a layer stack). No reference program: the checks are closed formulas, limits
// and two independently written algorithms (explicit summation of the
// multiply reflected beams, characteristic matrices).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <complex>
#include <memory>
#include <vector>

#include "ellipsometry/stack.hpp"
#include "optics/optics.hpp"

using namespace ellidft;
using namespace ellidft::ellipsometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const double degree = M_PI / 180.0;

std::shared_ptr<Material> constant(double n, double k = 0.0) { return std::make_shared<ConstantMaterial>(n, k); }

Stack bare(double n_ambient, Complex n_substrate) {
  return {constant(n_ambient), {}, constant(n_substrate.real(), n_substrate.imag())};
}

Complex q_of(Complex N, Complex kx) {
  Complex q = std::sqrt(N * N - kx * kx);
  return q.imag() < 0.0 ? -q : q;
}

// Second algorithm: characteristic matrices. In every medium the tangential
// field is U = a e^{i b} + b e^{-i b}, its partner V = eta (a e^{i b} - b e^{-i b})
// with eta = q for s polarisation (U: electric field) and eta = q / N^2 for
// p polarisation (U: magnetic field); q = N cos(theta). Across a layer of
// phase thickness beta = 2 pi d q / lambda
//   U(top) = cos(beta) U(bottom) - (i / eta) sin(beta) V(bottom)
//   V(top) = -i eta sin(beta) U(bottom) + cos(beta) V(bottom).
// The substrate carries only the transmitted wave, (U, V) = (1, eta_s); in the
// ambient U = 1 + r, V = eta_0 (1 - r).
Reflection matrix_method(const std::vector<Complex>& N, const std::vector<double>& thickness, double wavelength,
                         double angle_deg) {
  const Complex kx = N.front() * std::sin(angle_deg * degree);
  Reflection result;
  for (int polarisation = 0; polarisation < 2; ++polarisation) {
    auto eta = [&](std::size_t j) {
      const Complex q = q_of(N[j], kx);
      return polarisation == 0 ? q : q / (N[j] * N[j]);
    };
    Complex U = 1.0, V = eta(N.size() - 1);
    for (std::size_t j = N.size() - 2; j >= 1; --j) {
      const Complex beta = 2.0 * M_PI * thickness[j - 1] * q_of(N[j], kx) / wavelength;
      const Complex c = std::cos(beta), s = std::sin(beta), i(0.0, 1.0);
      const Complex U_top = c * U - i / eta(j) * s * V;
      const Complex V_top = -i * eta(j) * s * U + c * V;
      U = U_top;
      V = V_top;
    }
    const Complex r = (eta(0) * U - V) / (eta(0) * U + V);
    (polarisation == 0 ? result.rs : result.rp) = r;
  }
  return result;
}

Stack stack_of(const std::vector<Complex>& N, const std::vector<double>& thickness) {
  Stack stack;
  stack.ambient = constant(N.front().real(), N.front().imag());
  for (std::size_t j = 0; j < thickness.size(); ++j)
    stack.layers.push_back({constant(N[j + 1].real(), N[j + 1].imag()), thickness[j]});
  stack.substrate = constant(N.back().real(), N.back().imag());
  return stack;
}

}  // namespace

TEST_CASE("Fresnel reflection of a bare transparent substrate", "[ellipsometry]") {
  const double n1 = 1.0, n2 = 1.5;

  SECTION("normal incidence") {
    // R0 = |(n1 - n2) / (n1 + n2)|^2; r_p = -r_s in the convention used
    const Reflection r = reflect(bare(n1, n2), 600.0, 0.0);
    CHECK_THAT(std::norm(r.rs), WithinRel(0.04, 1e-14));
    CHECK_THAT(r.rs.real(), WithinRel(-0.2, 1e-14));
    CHECK_THAT(r.rp.real(), WithinRel(0.2, 1e-14));
    CHECK_THAT(psi_delta(r).psi_deg, WithinRel(45.0, 1e-14));
  }

  SECTION("oblique incidence, angles from Snell's law") {
    for (double angle : {20.0, 45.0, 70.0, 85.0}) {
      const double ci = std::cos(angle * degree);
      const double ct = std::cos(std::asin(n1 * std::sin(angle * degree) / n2));
      const double rs = (n1 * ci - n2 * ct) / (n1 * ci + n2 * ct);
      const double rp = (n2 * ci - n1 * ct) / (n2 * ci + n1 * ct);
      const Reflection r = reflect(bare(n1, n2), 600.0, angle);
      CHECK_THAT(r.rs.real(), WithinAbs(rs, 1e-15));
      CHECK_THAT(r.rp.real(), WithinAbs(rp, 1e-15));
      CHECK_THAT(r.rs.imag(), WithinAbs(0.0, 1e-15));
      CHECK_THAT(r.rp.imag(), WithinAbs(0.0, 1e-15));
    }
  }

  SECTION("Brewster angle") {
    // theta_B = arctan(n2 / n1): no reflection of p-polarised light
    const double brewster = std::atan(n2 / n1) / degree;
    const Reflection r = reflect(bare(n1, n2), 600.0, brewster);
    CHECK_THAT(std::abs(r.rp), WithinAbs(0.0, 1e-15));
    CHECK_THAT(psi_delta(r).psi_deg, WithinAbs(0.0, 1e-13));
    // Delta is 180 degrees below and 0 above the Brewster angle
    CHECK_THAT(psi_delta(bare(n1, n2), 600.0, brewster - 5.0).delta_deg, WithinAbs(180.0, 1e-12));
    const double above = psi_delta(bare(n1, n2), 600.0, brewster + 5.0).delta_deg;
    CHECK(std::min(above, 360.0 - above) < 1e-12);
  }

  SECTION("total internal reflection") {
    // from glass into air above the critical angle: |r| = 1 for both
    const Reflection r = reflect(bare(1.5, 1.0), 600.0, 60.0);
    CHECK_THAT(std::abs(r.rs), WithinRel(1.0, 1e-14));
    CHECK_THAT(std::abs(r.rp), WithinRel(1.0, 1e-14));
  }
}

TEST_CASE("an absorbing substrate is recovered from Psi and Delta", "[ellipsometry]") {
  // Eliminating cos(theta_t) from the two Fresnel formulas gives, for a bare
  // substrate in an ambient of index n0 and rho = r_p / r_s,
  //   N^2 = n0^2 sin^2(phi) [1 + tan^2(phi) ((1 - rho) / (1 + rho))^2].
  for (const Complex N : {Complex(3.9, 0.02), Complex(0.2, 3.4), Complex(1.5, 0.0), Complex(2.0, 1.0)}) {
    for (double angle : {50.0, 65.0, 75.0}) {
      const PsiDelta pd = psi_delta(bare(1.0, N), 500.0, angle);
      const Complex rho = std::tan(pd.psi_deg * degree) * std::exp(Complex(0.0, pd.delta_deg * degree));
      const double s = std::sin(angle * degree), t = std::tan(angle * degree);
      const Complex ratio = (1.0 - rho) / (1.0 + rho);
      const Complex eps = s * s * (1.0 + t * t * ratio * ratio);
      CHECK_THAT(std::abs(eps - N * N), WithinAbs(0.0, 1e-11 * std::norm(N)));
    }
  }
}

TEST_CASE("a single layer equals the sum of the multiply reflected beams", "[ellipsometry]") {
  // r = r01 + t01 t10 r12 e^{2 i beta} sum_k (r10 r12 e^{2 i beta})^k with the
  // Fresnel transmission coefficients
  //   t_s = 2 n1 cos_i / (n1 cos_i + n2 cos_t),  t_p = 2 n1 cos_i / (n2 cos_i + n1 cos_t)
  // (Wikipedia, "Fresnel equations").
  struct Case {
    Complex film, substrate;
    double thickness, wavelength, angle;
  };
  const Case cases[] = {{{1.59, 0.0}, {3.9, 0.02}, 120.0, 550.0, 70.0},
                        {{1.46, 0.0}, {3.9, 0.02}, 2.0, 633.0, 65.0},
                        {{1.7, 0.3}, {1.5, 0.0}, 300.0, 450.0, 55.0},
                        {{2.1, 0.0}, {1.45, 0.0}, 1000.0, 800.0, 30.0}};
  for (const Case& c : cases) {
    const Complex N0 = 1.0, kx = N0 * std::sin(c.angle * degree);
    const Complex q0 = q_of(N0, kx), q1 = q_of(c.film, kx), q2 = q_of(c.substrate, kx);
    auto rs = [](Complex qa, Complex qb) { return (qa - qb) / (qa + qb); };
    auto ts = [](Complex qa, Complex qb) { return 2.0 * qa / (qa + qb); };
    // with cos = q / N: r_p = (Nb^2 qa - Na^2 qb) / (...), t_p = 2 Na Nb qa / (Nb^2 qa + Na^2 qb)
    auto rp = [](Complex Na, Complex qa, Complex Nb, Complex qb) { return (Nb * Nb * qa - Na * Na * qb) / (Nb * Nb * qa + Na * Na * qb); };
    auto tp = [](Complex Na, Complex qa, Complex Nb, Complex qb) { return 2.0 * Na * Nb * qa / (Nb * Nb * qa + Na * Na * qb); };
    const Complex phase = std::exp(Complex(0.0, 2.0) * (2.0 * M_PI * c.thickness * q1 / c.wavelength));

    Complex sum_s = 0.0, sum_p = 0.0, term_s = 1.0, term_p = 1.0;
    for (int k = 0; k < 2000; ++k) {
      sum_s += term_s;
      sum_p += term_p;
      term_s *= rs(q1, q0) * rs(q1, q2) * phase;
      term_p *= rp(c.film, q1, N0, q0) * rp(c.film, q1, c.substrate, q2) * phase;
    }
    const Complex series_s = rs(q0, q1) + ts(q0, q1) * ts(q1, q0) * rs(q1, q2) * phase * sum_s;
    const Complex series_p = rp(N0, q0, c.film, q1) + tp(N0, q0, c.film, q1) * tp(c.film, q1, N0, q0) * rp(c.film, q1, c.substrate, q2) * phase * sum_p;

    const Reflection r = reflect(stack_of({N0, c.film, c.substrate}, {c.thickness}), c.wavelength, c.angle);
    CHECK_THAT(std::abs(r.rs - series_s), WithinAbs(0.0, 1e-13));
    CHECK_THAT(std::abs(r.rp - series_p), WithinAbs(0.0, 1e-13));
  }
}

TEST_CASE("a glass plate in air at normal incidence", "[ellipsometry]") {
  // Wikipedia, "Transfer-matrix method (optics)": for a single layer of index
  // n and thickness d in air, with k' = 2 pi n / lambda,
  //   r = (1/n - n) sin(k'd) / ((n + 1/n) sin(k'd) + 2 i cos(k'd)).
  // Only the modulus is compared (the page does not state its phase convention).
  const double n = 1.5, wavelength = 600.0;
  for (double d : {0.0, 37.0, 100.0, 250.0, 1234.5}) {
    const double kd = 2.0 * M_PI * n * d / wavelength;
    const Complex r_page = (1.0 / n - n) * std::sin(kd) / Complex((n + 1.0 / n) * std::sin(kd), 2.0 * std::cos(kd));
    const Reflection r = reflect(stack_of({1.0, n, 1.0}, {d}), wavelength, 0.0);
    CHECK_THAT(std::abs(r.rs), WithinAbs(std::abs(r_page), 1e-14));
    CHECK_THAT(std::abs(r.rp), WithinAbs(std::abs(r_page), 1e-14));
  }
  // no reflection when k'd is a multiple of pi
  const Reflection r = reflect(stack_of({1.0, n, 1.0}, {wavelength / (2.0 * n)}), wavelength, 0.0);
  CHECK_THAT(std::abs(r.rs), WithinAbs(0.0, 1e-15));
}

TEST_CASE("layer stacks: limits and consistency", "[ellipsometry]") {
  const Complex substrate(3.9, 0.02);
  const Reflection r_bare = reflect(bare(1.0, substrate), 550.0, 70.0);

  SECTION("a layer of zero thickness has no effect") {
    const Reflection r = reflect(stack_of({1.0, {1.5, 0.1}, substrate}, {0.0}), 550.0, 70.0);
    CHECK_THAT(std::abs(r.rs - r_bare.rs), WithinAbs(0.0, 1e-15));
    CHECK_THAT(std::abs(r.rp - r_bare.rp), WithinAbs(0.0, 1e-15));
  }
  SECTION("a layer of the substrate material has no effect") {
    const Reflection r = reflect(stack_of({1.0, substrate, substrate}, {80.0}), 550.0, 70.0);
    CHECK_THAT(std::abs(r.rs - r_bare.rs), WithinAbs(0.0, 1e-15));
    CHECK_THAT(std::abs(r.rp - r_bare.rp), WithinAbs(0.0, 1e-15));
  }
  SECTION("a layer of the ambient material only shifts the common phase") {
    const Stack stack = stack_of({1.0, 1.0, substrate}, {300.0});
    const PsiDelta with = psi_delta(stack, 550.0, 70.0), without = psi_delta(r_bare);
    CHECK_THAT(with.psi_deg, WithinAbs(without.psi_deg, 1e-12));
    CHECK_THAT(with.delta_deg, WithinAbs(without.delta_deg, 1e-12));
  }
  SECTION("a transparent layer with phase thickness pi is absent") {
    const double n = 1.6, wavelength = 550.0, angle = 70.0;
    const double q = std::sqrt(n * n - std::sin(angle * degree) * std::sin(angle * degree));
    const Reflection r = reflect(stack_of({1.0, n, substrate}, {wavelength / (2.0 * q)}), wavelength, angle);
    CHECK_THAT(std::abs(r.rs - r_bare.rs), WithinAbs(0.0, 1e-14));
    CHECK_THAT(std::abs(r.rp - r_bare.rp), WithinAbs(0.0, 1e-14));
  }
  SECTION("splitting a layer in two changes nothing") {
    const Complex film(1.7, 0.05);
    const Reflection one = reflect(stack_of({1.0, film, substrate}, {137.0}), 480.0, 62.0);
    const Reflection two = reflect(stack_of({1.0, film, film, substrate}, {50.0, 87.0}), 480.0, 62.0);
    CHECK_THAT(std::abs(one.rs - two.rs), WithinAbs(0.0, 1e-15));
    CHECK_THAT(std::abs(one.rp - two.rp), WithinAbs(0.0, 1e-15));
  }
  SECTION("several layers: recursion against characteristic matrices") {
    const std::vector<Complex> N = {1.0, {1.38, 0.0}, {2.05, 0.01}, {1.59, 0.0}, {1.46, 0.0}, {3.9, 0.02}};
    const std::vector<double> d = {90.0, 45.0, 310.0, 2.5};
    double largest = 0.0;
    for (double wavelength : {350.0, 500.0, 633.0, 900.0}) {
      for (double angle : {0.0, 45.0, 65.0, 75.0}) {
        const Reflection a = reflect(stack_of(N, d), wavelength, angle), b = matrix_method(N, d, wavelength, angle);
        largest = std::max({largest, std::abs(a.rs - b.rs), std::abs(a.rp - b.rp)});
      }
    }
    CHECK(largest < 1e-13);
    SUCCEED("largest difference between the two algorithms: " << largest);
  }
  SECTION("a lossless stack on a transparent substrate reflects at most everything") {
    const std::vector<Complex> N = {1.0, 2.3, 1.38, 2.3, 1.38, 1.52};
    const std::vector<double> d = {60.0, 100.0, 60.0, 100.0};
    for (double wavelength = 300.0; wavelength <= 900.0; wavelength += 25.0) {
      const Reflection r = reflect(stack_of(N, d), wavelength, 50.0);
      CHECK(std::norm(r.rs) < 1.0);
      CHECK(std::norm(r.rp) < 1.0);
    }
  }
}

TEST_CASE("forward model rejects invalid input", "[ellipsometry]") {
  Stack stack = bare(1.0, 1.5);
  CHECK_THROWS_AS(reflect(stack, 0.0, 60.0), std::invalid_argument);
  CHECK_THROWS_AS(reflect(stack, 500.0, 90.0), std::invalid_argument);
  CHECK_THROWS_AS(reflect(stack, 500.0, -1.0), std::invalid_argument);
  stack.layers.push_back({constant(1.5), -1.0});
  CHECK_THROWS_AS(reflect(stack, 500.0, 60.0), std::invalid_argument);
  stack.layers.back() = {nullptr, 10.0};
  CHECK_THROWS_AS(reflect(stack, 500.0, 60.0), std::invalid_argument);
  CHECK_THROWS_AS(reflect(Stack{}, 500.0, 60.0), std::invalid_argument);
}

TEST_CASE("material models", "[ellipsometry]") {
  SECTION("Cauchy") {
    CauchyMaterial m(1.5, 0.006, 0.0002);
    // lambda = 0.5 um: 1.5 + 0.006 / 0.25 + 0.0002 / 0.0625
    CHECK_THAT(m.index(500.0).real(), WithinRel(1.5 + 0.024 + 0.0032, 1e-15));
    CHECK(m.index(500.0).imag() == 0.0);
    CHECK(m.parameter_names() == std::vector<std::string>{"A", "B", "C"});
    m.set_parameter(1, 0.0);
    m.set_parameter(2, 0.0);
    CHECK(m.index(321.0).real() == 1.5);
  }
  SECTION("Sellmeier") {
    // Polystyrene, one term: B = 1.4435, C = 0.020216 um^2. The value at
    // 589 nm is that of tests/reference/data/m6/guide_values.json.
    SellmeierMaterial ps({1.4435}, {0.020216});
    CHECK_THAT(ps.index(589.0).real(), WithinRel(1.591484070596161, 1e-14));
    CHECK(ps.parameter_names() == std::vector<std::string>{"B1", "C1"});
    CHECK(ps.parameter(0) == 1.4435);
    CHECK(ps.parameter(1) == 0.020216);
    // two terms add up
    SellmeierMaterial two({0.7, 0.5}, {0.01, 0.02});
    const double l2 = 0.36;
    CHECK_THAT(two.index(600.0).real(), WithinRel(std::sqrt(1.0 + 0.7 * l2 / (l2 - 0.01) + 0.5 * l2 / (l2 - 0.02)), 1e-15));
    // inside the absorption (lambda^2 just below C) n^2 is negative
    CHECK_THROWS_AS(ps.index(140.0), std::domain_error);
    CHECK_THROWS_AS(SellmeierMaterial({1.0}, {}), std::invalid_argument);
  }
  SECTION("table") {
    TabulatedMaterial m({400.0, 500.0, 700.0}, {{1.6, 0.2}, {1.5, 0.1}, {1.4, 0.0}});
    CHECK(m.index(400.0) == Complex(1.6, 0.2));
    CHECK(m.index(700.0) == Complex(1.4, 0.0));
    CHECK_THAT(m.index(450.0).real(), WithinRel(1.55, 1e-15));
    CHECK_THAT(m.index(450.0).imag(), WithinRel(0.15, 1e-15));
    CHECK_THAT(m.index(650.0).real(), WithinRel(1.425, 1e-15));
    CHECK_THROWS_AS(m.index(399.9), std::out_of_range);   // no extrapolation
    CHECK_THROWS_AS(m.index(700.1), std::out_of_range);
    CHECK_THROWS_AS(TabulatedMaterial({500.0, 400.0}, {{1.0, 0.0}, {1.0, 0.0}}), std::invalid_argument);
    CHECK_THROWS_AS(TabulatedMaterial({400.0}, {}), std::invalid_argument);
  }
  SECTION("polarisability and density (Lorentz-Lorenz)") {
    const std::vector<double> wavelengths = {400.0, 600.0, 800.0};
    const std::vector<Complex> alpha = {{14.0, 0.3}, {13.3, 0.0}, {13.1, 0.0}};
    LorentzLorenzMaterial m(wavelengths, alpha, 104.15, 1.05);
    const Complex expected = optics::lorentz_lorenz(Complex(13.3, 0.0), optics::number_density(1.05, 104.15));
    CHECK_THAT(std::abs(m.index(600.0) - expected), WithinAbs(0.0, 1e-15));
    CHECK(m.index(400.0).imag() > 0.0);
    // the density is the adjustable parameter; a denser material refracts more
    CHECK(m.parameter_names() == std::vector<std::string>{"density"});
    const double before = m.index(600.0).real();
    m.set_parameter(0, 1.10);
    CHECK(m.parameter(0) == 1.10);
    CHECK(m.index(600.0).real() > before);
    // the inverse relation returns the polarisability
    const Complex back = optics::lorentz_lorenz_inverse(m.index(600.0), optics::number_density(1.10, 104.15));
    CHECK_THAT(std::abs(back - Complex(13.3, 0.0)), WithinAbs(0.0, 1e-13));
  }
  SECTION("materials without parameters") {
    TabulatedMaterial m({400.0, 500.0}, {{1.0, 0.0}, {1.0, 0.0}});
    CHECK(m.parameter_names().empty());
    CHECK_THROWS_AS(m.parameter(0), std::out_of_range);
    CHECK_THROWS_AS(m.set_parameter(0, 1.0), std::out_of_range);
  }
}
