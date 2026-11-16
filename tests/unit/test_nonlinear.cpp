// Frequency-dependent non-linear response from real-time propagation.
// Independent routes used for comparison:
//   - linear response in a static field (SCF with the field, then the
//     response equations): d alpha(w) / dF = beta(-w; w, 0),
//     d^2 alpha(w) / dF^2 = gamma(-w; w, 0, 0);
//   - the low-frequency expansion: for all fields and the response along one
//     direction, permutation symmetry leaves one quadratic term,
//     X(w_sigma; w_1, ...) = X(0) + A (w_sigma^2 + w_1^2 + ...) + O(w^4),
//     so the dispersion of second-harmonic generation is 3 times that of the
//     Pockels effect, that of the optical Kerr effect 2 times and that of
//     third-harmonic generation 6 times that of the dc Kerr effect;
//   - excitation energies of states without oscillator strength from the
//     response code: two-photon absorption must peak at half of them.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdio>
#include <string>

#include "core/xyz.hpp"
#include "optics/optics.hpp"
#include "realtime/nonlinear.hpp"
#include "response/excitations.hpp"
#include "response/polarizability.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Molecule h2() { return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2.xyz"); }

Molecule tilted_water() {
  Molecule mol = read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2o.xyz");
  const Eigen::Vector3d u = Eigen::Vector3d(0.4, -0.9, 0.6).normalized();
  Eigen::Matrix3d cross;
  cross << 0.0, -u.z(), u.y(), u.z(), 0.0, -u.x(), -u.y(), u.x(), 0.0;
  const Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + std::sin(0.77) * cross + (1.0 - std::cos(0.77)) * cross * cross;
  for (Atom& atom : mol.atoms) {
    const Eigen::Vector3d r = R * Eigen::Vector3d(atom.r[0], atom.r[1], atom.r[2]);
    atom.r = {r.x(), r.y(), r.z()};
  }
  return mol;
}

ScfOptions tight() {
  ScfOptions options;
  options.energy_tolerance = 1e-12;
  options.error_tolerance = 1e-10;
  return options;
}

double dot(const Vec3& a, const Eigen::Vector3d& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// u.alpha(w).u from the response equations for the molecule in the static field F u.
double alpha_in_field(const Molecule& molecule, const BasisSet& basis, const Method& method, JKBuilder& jk,
                      const Eigen::Vector3d& u, double F, double omega) {
  ScfOptions options = tight();
  options.electric_field = {F * u[0], F * u[1], F * u[2]};
  const ScfResult scf = run_scf(molecule, basis, method, jk, options);
  REQUIRE(scf.converged);
  const ResponseOperator op(molecule, basis, method, scf, jk);
  LinearResponseOptions response;
  response.tolerance = 1e-11;
  const Polarizability alpha = polarizability(op, omega, 0.0, response);
  REQUIRE(alpha.converged);
  return u.dot(alpha.real * u);
}

// u.mu for the molecule in the static field F u.
double dipole_in_field(const Molecule& molecule, const BasisSet& basis, const Method& method, JKBuilder& jk,
                       const Eigen::Vector3d& u, double F) {
  ScfOptions options = tight();
  options.electric_field = {F * u[0], F * u[1], F * u[2]};
  const ScfResult scf = run_scf(molecule, basis, method, jk, options);
  REQUIRE(scf.converged);
  return dot(scf.dipole, u);
}

}  // namespace

namespace {

void check_second_order_of_water(const Method& method) {
  INFO(method.functional.name);
  const Molecule molecule = tilted_water();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const Eigen::Vector3d u = Eigen::Vector3d(1.0, -2.0, 1.5).normalized();
  const double w = 0.08;

  // Independent values: static beta from dipoles, Pockels from alpha(w), both in static fields.
  const double h = 0.004;
  const double beta_static = (dipole_in_field(molecule, basis, method, jk, u, h) - 2.0 * dot(scf.dipole, u) +
                              dipole_in_field(molecule, basis, method, jk, u, -h)) / (h * h);
  const double alpha_w = alpha_in_field(molecule, basis, method, jk, u, 0.0, w);
  const double pockels_reference = (alpha_in_field(molecule, basis, method, jk, u, h, w) -
                                    alpha_in_field(molecule, basis, method, jk, u, -h, w)) / (2.0 * h);

  HarmonicFieldOptions options;
  options.omega = w;
  options.optical = {u[0], u[1], u[2]};
  options.static_field = {u[0], u[1], u[2]};
  options.amplitude = 0.01;
  options.ramp_time = 160.0;
  options.n_periods = 1;
  options.propagation.time_step = 0.2;
  options.propagation.fock_tolerance = 1e-11;
  const HarmonicResponse r = harmonic_response(propagator, options);
  REQUIRE(r.steps_converged);

  const double alpha_rt = dot(r.cosine[0][1], u);
  const double shg = 4.0 * dot(r.cosine[1][2], u);        // beta(-2w; w, w)
  const double pockels = dot(r.cosine[1][1], u);          // beta(-w; w, 0)
  std::printf("non-linear, water, %s, w = %.2f: alpha %.6f (response code %.6f); beta static %.5f, Pockels %.5f "
              "(alpha in static fields %.5f), second harmonic %.5f; fit residuals %.1e %.1e; %d Fock builds\n",
              method.functional.name.c_str(), w, alpha_rt, alpha_w, beta_static, pockels, pockels_reference, shg, r.residual[0], r.residual[1], r.fock_builds);
  CHECK_THAT(alpha_rt, WithinRel(alpha_w, 1e-4));
  CHECK_THAT(pockels, WithinRel(pockels_reference, 2e-3));
  // dispersion: second harmonic 3 times the Pockels effect, up to O(w^4)
  const double ratio = (shg - beta_static) / (pockels_reference - beta_static);
  std::printf("non-linear, water: dispersion ratio second harmonic / Pockels %.3f (3 for w -> 0)\n", ratio);
  CHECK(std::abs(pockels_reference - beta_static) > 1e-3 * std::abs(beta_static));   // there is dispersion to compare
  CHECK_THAT(ratio, WithinAbs(3.0, 0.3));
  // no absorption below the first excitation: nothing in quadrature
  CHECK(std::abs(dot(r.sine[0][1], u)) < 1e-3 * std::abs(alpha_rt));
  CHECK(r.residual[0] < 1e-3);
}

}  // namespace

TEST_CASE("non-linear: second order of water (second harmonic, Pockels effect)", "[nonlinear][slow]") {
  check_second_order_of_water(Method{});
}

TEST_CASE("non-linear: second order of water with density functionals", "[nonlinear][slow]") {
  // small grid: the propagation and the response code use the same one
  for (const char* functional : {"PBE0", "CAM-B3LYP"}) {
    Method method;
    method.functional = Functional::by_name(functional);
    method.grid = {20, 50};
    check_second_order_of_water(method);
  }
}

TEST_CASE("non-linear: third order of H2 (optical and dc Kerr effect, third harmonic)", "[nonlinear][slow]") {
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("cc-pVDZ"));
  Method method;
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const Eigen::Vector3d z(0.0, 0.0, 1.0);
  const double w = 0.1;

  // Independent values from calculations in static fields.
  const double h = 0.005;
  const double m1 = dipole_in_field(molecule, basis, method, jk, z, h), m2 = dipole_in_field(molecule, basis, method, jk, z, 2.0 * h);
  const double gamma_static = (m2 - 2.0 * m1) / (h * h * h);   // odd function: mu = alpha F + gamma F^3 / 6
  const double a0 = alpha_in_field(molecule, basis, method, jk, z, 0.0, w), a1 = alpha_in_field(molecule, basis, method, jk, z, h, w),
               a2 = alpha_in_field(molecule, basis, method, jk, z, 2.0 * h, w);
  // alpha(F) = a0 + gamma F^2 / 2 + O(F^4): eliminate the F^4 term
  const double dc_kerr_reference = 2.0 * (16.0 * (a1 - a0) - (a2 - a0)) / (12.0 * h * h);

  HarmonicFieldOptions options;
  options.omega = w;
  options.optical = {0.0, 0.0, 1.0};
  options.amplitude = 0.01;      // the fifth order enters the third as 5 a^2 mu5
  options.ramp_time = 300.0;
  options.n_periods = 2;
  options.propagation.time_step = 0.2;
  options.propagation.fock_tolerance = 1e-12;
  const HarmonicResponse optical = harmonic_response(propagator, options);
  const double f0 = 0.5;         // static part: half the optical amplitude
  options.static_field = {0.0, 0.0, f0};
  const HarmonicResponse with_static = harmonic_response(propagator, options);
  REQUIRE(optical.steps_converged);
  REQUIRE(with_static.steps_converged);

  const double kerr = 8.0 * optical.cosine[2][1][2];             // gamma(-w; w, w, -w)
  const double thg = 24.0 * optical.cosine[2][3][2];             // gamma(-3w; w, w, w)
  const double dc_kerr = 2.0 * (with_static.cosine[2][1][2] - optical.cosine[2][1][2]) / (f0 * f0);   // gamma(-w; w, 0, 0)
  std::printf("non-linear, H2, w = %.2f: gamma static %.3f, dc Kerr %.3f (alpha in static fields %.3f), optical Kerr %.3f, "
              "third harmonic %.3f; fit residual third order %.1e\n",
              w, gamma_static, dc_kerr, dc_kerr_reference, kerr, thg, optical.residual[2]);
  std::printf("non-linear, H2: quadrature part of the optical Kerr term %.3f, %d + %d Fock builds\n",
              8.0 * optical.sine[2][1][2], optical.fock_builds, with_static.fock_builds);
  CHECK_THAT(dc_kerr, WithinRel(dc_kerr_reference, 1e-2));
  const double kerr_ratio = (kerr - gamma_static) / (dc_kerr_reference - gamma_static);
  const double thg_ratio = (thg - gamma_static) / (dc_kerr_reference - gamma_static);
  std::printf("non-linear, H2: dispersion ratios optical Kerr / dc Kerr %.3f (2 for w -> 0), third harmonic / dc Kerr %.3f (6)\n",
              kerr_ratio, thg_ratio);
  CHECK_THAT(kerr_ratio, WithinAbs(2.0, 0.3));
  CHECK_THAT(thg_ratio, WithinAbs(6.0, 1.5));
  // centre of symmetry: no second order
  CHECK(std::abs(optical.cosine[1][2][2]) < 1e-6);
}

TEST_CASE("non-linear: two-photon absorption of H2 peaks at half the energy of a dark state", "[nonlinear][slow]") {
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("cc-pVDZ"));
  Method method;
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const ResponseOperator op(molecule, basis, method, scf, jk);
  ExcitationOptions excitation_options;
  excitation_options.n_states = 8;
  const ExcitationResult states = excitations(op, excitation_options);
  REQUIRE(states.converged);
  // lowest allowed state along z, and the dark states below 1.2 Hartree
  double first_allowed = 0.0;
  for (const Excitation& s : states.states)
    if (first_allowed == 0.0 && std::abs(s.transition_dipole[2]) > 1e-3) first_allowed = s.energy;
  REQUIRE(first_allowed > 0.0);
  std::printf("non-linear, H2: excitations");
  for (const Excitation& s : states.states) std::printf("  %.4f (f = %.3f)", s.energy, s.oscillator_strength);
  std::printf("\n");

  TwoPhotonOptions options;
  options.amplitude = 0.01;
  options.duration = 300.0;
  options.propagation.propagator = Propagator::Composition4;
  options.propagation.time_step = 0.2;
  options.propagation.fock_tolerance = 1e-11;
  options.propagation.max_fock_builds_per_step = 20;
  int n_tested = 0;
  for (const Excitation& s : states.states) {
    const bool dark = s.oscillator_strength < 1e-8;
    if (!dark || 0.5 * s.energy >= first_allowed - 0.1) continue;
    double on = 0.0, off = 0.0;
    for (double shift : {0.0, 0.06}) {
      options.omega = 0.5 * s.energy - shift;
      const TwoPhotonResult r = two_photon_absorption(propagator, options);
      REQUIRE(r.steps_converged);
      (shift == 0.0 ? on : off) = r.cross_section_gm;
      std::printf("non-linear, H2, dark state %.4f: w = %.4f  W2 = %.3e  W4 = %.3e  sigma2 = %.4e a.u. = %.4e GM\n", s.energy,
                  options.omega, r.one_photon_coefficient, r.two_photon_coefficient, r.cross_section_au, r.cross_section_gm);
    }
    // strong at half the energy of the dark state, nothing 0.06 Hartree below
    CHECK(on > 0.0);
    CHECK(std::abs(off) < 1e-3 * on);
    ++n_tested;
  }
  CHECK(n_tested >= 1);
}

// Short version for every test run: H2 in the minimal basis (two orbitals).
TEST_CASE("non-linear: harmonic analysis and pulse absorption on a minimal example", "[nonlinear]") {
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  Method method;
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const Eigen::Vector3d z(0.0, 0.0, 1.0);
  const double w = 0.15;

  HarmonicFieldOptions options;
  options.omega = w;
  options.optical = {0.0, 0.0, 1.0};
  options.static_field = {0.0, 0.0, 0.5};
  options.amplitude = 0.01;
  options.ramp_time = 150.0;
  options.n_periods = 1;
  options.propagation.time_step = 0.2;
  const HarmonicResponse r = harmonic_response(propagator, options);
  REQUIRE(r.steps_converged);
  // first order: alpha(w) at the optical frequency, alpha(0) times the static part
  CHECK_THAT(r.cosine[0][1][2], WithinRel(alpha_in_field(molecule, basis, method, jk, z, 0.0, w), 1e-4));
  CHECK_THAT(r.cosine[0][0][2], WithinRel(0.5 * alpha_in_field(molecule, basis, method, jk, z, 0.0, 0.0), 1e-4));
  CHECK(std::abs(r.sine[0][1][2]) < 1e-4 * std::abs(r.cosine[0][1][2]));   // no absorption
  // centre of symmetry: no second order; nothing perpendicular to the bond
  for (std::size_t m = 0; m < 4; ++m) {
    CHECK(std::abs(r.cosine[1][m][2]) < 1e-6);
    CHECK(std::abs(r.cosine[0][m][0]) < 1e-12);
  }
  CHECK(r.residual[0] < 1e-3);

  // a pulse far below the only excitation leaves next to nothing
  TwoPhotonOptions pulse;
  pulse.omega = 0.2;
  pulse.duration = 100.0;
  pulse.propagation.propagator = Propagator::Composition4;
  pulse.propagation.time_step = 0.2;
  const TwoPhotonResult tp = two_photon_absorption(propagator, pulse);
  CHECK(tp.steps_converged);
  CHECK(std::abs(tp.one_photon_coefficient) < 1e-2);
  // conversion: 1 a.u. of the cross-section is Bohr^4 times the atomic unit of time
  if (tp.cross_section_au != 0.0)
    CHECK_THAT(tp.cross_section_gm / tp.cross_section_au,
               WithinRel(std::pow(0.529177210903e-8, 4) * 2.4188843265864e-17 / 1e-50, 1e-12));
}

TEST_CASE("non-linear: invalid requests", "[nonlinear]") {
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  Method method;
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  HarmonicFieldOptions harmonic;
  CHECK_THROWS_AS(harmonic_response(propagator, harmonic), std::invalid_argument);   // omega missing
  harmonic.omega = 0.1;
  harmonic.n_periods = 0;
  CHECK_THROWS_AS(harmonic_response(propagator, harmonic), std::invalid_argument);
  TwoPhotonOptions two_photon;
  CHECK_THROWS_AS(two_photon_absorption(propagator, two_photon), std::invalid_argument);   // omega missing
  two_photon.omega = 0.1;
  two_photon.polarization = {0.0, 0.0, 0.0};
  CHECK_THROWS_AS(two_photon_absorption(propagator, two_photon), std::invalid_argument);
  CHECK(smooth_switch(-1.0) == 0.0);
  CHECK(smooth_switch(2.0) == 1.0);
  CHECK_THAT(smooth_switch(0.5), WithinRel(0.5, 1e-14));
}

// ---------------------------------------------------------------------------
// Quantities for disordered materials
// ---------------------------------------------------------------------------

TEST_CASE("non-linear: seven directions average polynomials of degree four exactly", "[nonlinear]") {
  const auto directions = averaging_directions();
  auto average = [&](int px, int py, int pz) {
    double sum = 0.0;
    for (const AveragingDirection& d : directions)
      sum += d.weight * std::pow(d.u[0], px) * std::pow(d.u[1], py) * std::pow(d.u[2], pz);
    return sum;
  };
  // averages over the unit sphere: <1> = 1, <x^2> = 1/3, <x^4> = 1/5, <x^2 y^2> = 1/15, odd powers 0
  CHECK_THAT(average(0, 0, 0), WithinRel(1.0, 1e-14));
  for (int k = 0; k < 3; ++k) {
    CHECK_THAT(average(k == 0 ? 2 : 0, k == 1 ? 2 : 0, k == 2 ? 2 : 0), WithinRel(1.0 / 3.0, 1e-14));
    CHECK_THAT(average(k == 0 ? 4 : 0, k == 1 ? 4 : 0, k == 2 ? 4 : 0), WithinRel(0.2, 1e-14));
    CHECK_THAT(average(k == 0 ? 0 : 2, k == 1 ? 0 : 2, k == 2 ? 0 : 2), WithinRel(1.0 / 15.0, 1e-14));
  }
  CHECK_THAT(average(1, 1, 0), WithinAbs(0.0, 1e-15));
  CHECK_THAT(average(3, 1, 0), WithinAbs(0.0, 1e-15));
  CHECK_THAT(average(2, 1, 1), WithinAbs(0.0, 1e-15));
  CHECK_THAT(average(1, 1, 2), WithinAbs(0.0, 1e-15));
  // a tensor of rank four: sum of w u.g.uuu = (1/15) sum_ij (g_iijj + g_ijij + g_ijji)
  double g[3][3][3][3];
  unsigned state = 12345u;
  for (auto& a : g) for (auto& b : a) for (auto& c : b) for (double& x : c) {
    state = state * 1664525u + 1013904223u;
    x = double(state >> 8) / double(1u << 24) - 0.5;
  }
  double from_directions = 0.0, closed = 0.0;
  for (const AveragingDirection& d : directions)
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l)
      from_directions += d.weight * g[i][j][k][l] * d.u[static_cast<std::size_t>(i)] * d.u[static_cast<std::size_t>(j)] *
                         d.u[static_cast<std::size_t>(k)] * d.u[static_cast<std::size_t>(l)];
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) closed += (g[i][i][j][j] + g[i][j][i][j] + g[i][j][j][i]) / 15.0;
  CHECK_THAT(from_directions, WithinAbs(closed, 1e-14));
}

TEST_CASE("non-linear: the index n2 agrees with the Lorentz-Lorenz relation for a field-dependent polarisability",
          "[nonlinear]") {
  // Material: polarisability volume 13.3 A^3 per unit, density of units as in polystyrene; gamma invented.
  const double alpha_A3 = 13.3, N = optics::number_density(1.05, 104.15), gamma = 2.5e4;
  const double n0 = optics::lorentz_lorenz(alpha_A3, N).real();
  const NonlinearIndex index = nonlinear_index(gamma, N, n0);
  CHECK_THAT(index.local_field_factor, WithinRel((n0 * n0 + 2.0) / 3.0, 1e-15));
  // Field amplitude E in the material, local field L E; the polarisability at the optical frequency
  // becomes alpha + gamma (L E)^2 / 8 (atomic units; 1 a.u. = 0.148184711 A^3).
  const double E = 1e-4, L = index.local_field_factor;
  const double alpha_shifted = alpha_A3 + optics::angstrom3_per_au * gamma * (L * E) * (L * E) / 8.0;
  const double delta_n = optics::lorentz_lorenz(alpha_shifted, N).real() - n0;
  const double intensity = n0 * 137.036 * E * E / (8.0 * M_PI);   // atomic units
  CHECK_THAT(index.n2_au * intensity, WithinRel(delta_n, 1e-5));
  CHECK_THAT(index.n2_cm2_per_W, WithinRel(1e4 * index.n2_m2_per_W, 1e-14));
  CHECK(index.n2_m2_per_W > 0.0);
  CHECK_THROWS_AS(nonlinear_index(gamma, 0.0, 1.5), std::invalid_argument);
  CHECK_THROWS_AS(nonlinear_index(gamma, N, 0.5), std::invalid_argument);
}

TEST_CASE("non-linear: transforms of the pulse", "[nonlinear]") {
  const double T = 240.0;
  CHECK_THAT(envelope_squared_transform(T, 0.0).real(), WithinRel(3.0 * T / 8.0, 1e-10));   // integral of sin^4
  CHECK_THAT(envelope_squared_transform(T, 0.0).imag(), WithinAbs(0.0, 1e-10));
  CHECK(std::abs(envelope_squared_transform(T, 0.5)) < 1e-3 * T);                           // far off resonance
  // on resonance, many cycles: half of the integral of sin^2, i.e. T / 4
  CHECK_THAT(std::abs(pulse_transform(0.4, T, 0.4)), WithinRel(T / 4.0, 1e-3));
  CHECK_THAT(std::abs(pulse_transform(0.4, T, -0.4)), WithinRel(T / 4.0, 1e-3));
  CHECK(std::abs(pulse_transform(0.4, T, 0.9)) < 1e-3 * T);
}

TEST_CASE("non-linear: the energy left by a pulse equals the oscillator strengths times the pulse spectrum", "[nonlinear]") {
  // First order: W = F0^2 sum_f w_f |e.t_f|^2 |integral f(t) e^{i w_f t} dt|^2, with the excitation energies and
  // transition dipoles of the response code. Checks pulse, units and the separation of the orders.
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  Method method;
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const ResponseOperator op(molecule, basis, method, scf, jk);
  ExcitationOptions excitation_options;
  excitation_options.n_states = 1;
  const ExcitationResult states = excitations(op, excitation_options);
  REQUIRE(states.converged);
  const Excitation& s = states.states.front();

  TwoPhotonOptions pulse;
  pulse.omega = s.energy - 0.03;   // slightly off the line
  pulse.duration = 100.0;
  pulse.amplitude = 0.002;
  pulse.propagation.propagator = Propagator::Composition4;
  pulse.propagation.time_step = 0.2;
  pulse.propagation.fock_tolerance = 1e-12;
  const TwoPhotonResult r = two_photon_absorption(propagator, pulse);
  REQUIRE(r.steps_converged);
  const double expected = s.energy * s.transition_dipole[2] * s.transition_dipole[2] *
                          std::norm(pulse_transform(pulse.omega, pulse.duration, s.energy));
  std::printf("non-linear, H2 minimal basis: one-photon absorption W2 = %.6f, from the response code %.6f\n",
              r.one_photon_coefficient, expected);
  CHECK_THAT(r.one_photon_coefficient, WithinRel(expected, 1e-3));
}

TEST_CASE("non-linear: absorption by a dark state in a static field", "[nonlinear][slow]") {
  // A static field mixes a dark state with allowed ones. The response code for the molecule in the field gives
  // the induced strength; the pulse absorption of the propagation in the same field must agree. This is the
  // limit of two-photon absorption in which one of the two photons has zero frequency.
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("cc-pVDZ"));
  Method method;
  DirectJK jk(basis);
  const double Fs = 0.02;
  ScfOptions in_field = tight();
  in_field.electric_field = {0.0, 0.0, Fs};
  const ScfResult scf = run_scf(molecule, basis, method, jk, in_field);
  REQUIRE(scf.converged);
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const ResponseOperator op(molecule, basis, method, scf, jk);
  ExcitationOptions excitation_options;
  excitation_options.n_states = 9;   // all of them: 1 occupied, 9 virtual orbitals
  const ExcitationResult states = excitations(op, excitation_options);
  REQUIRE(states.converged);
  const Excitation& dark = states.states[1];   // 0.78 Hartree without field, no oscillator strength there

  TwoPhotonOptions pulse;
  pulse.omega = dark.energy;
  pulse.duration = 200.0;
  pulse.amplitude = 0.002;
  pulse.static_field = {0.0, 0.0, Fs};
  pulse.propagation.propagator = Propagator::Composition4;
  pulse.propagation.time_step = 0.2;
  pulse.propagation.fock_tolerance = 1e-12;
  pulse.propagation.max_fock_builds_per_step = 20;
  const TwoPhotonResult r = two_photon_absorption(propagator, pulse);
  REQUIRE(r.steps_converged);
  double expected = 0.0, from_dark_state = 0.0;
  for (const Excitation& s : states.states) {
    const double part = s.energy * s.transition_dipole[2] * s.transition_dipole[2] *
                        std::norm(pulse_transform(pulse.omega, pulse.duration, s.energy));
    expected += part;
    if (&s == &dark) from_dark_state = part;
  }
  std::printf("non-linear, H2 in a static field %.2f: dark state at %.4f with induced transition dipole %.5f; absorption "
              "W2 = %.6f, from the response code %.6f (of which the dark state %.6f)\n",
              Fs, dark.energy, dark.transition_dipole[2], r.one_photon_coefficient, expected, from_dark_state);
  CHECK(from_dark_state > 0.9 * expected);   // the induced line dominates at this frequency
  CHECK_THAT(r.one_photon_coefficient, WithinRel(expected, 1e-2));
}

TEST_CASE("non-linear: orientational averages do not depend on the orientation of the molecule", "[nonlinear][slow]") {
  // H2 with s functions only responds along the bond alone: for any orientation of the molecule the average
  // over orientations is one fifth (<cos^4> = 1/5) of the value along the bond.
  const Eigen::Vector3d axis = Eigen::Vector3d(0.3, -0.5, 0.8).normalized();
  Molecule aligned = h2(), tilted = h2();
  const double half = 0.5 * std::abs(aligned.atoms[1].r[2] - aligned.atoms[0].r[2]);
  for (int k = 0; k < 2; ++k)
    for (std::size_t c = 0; c < 3; ++c) tilted.atoms[static_cast<std::size_t>(k)].r[c] = (k == 0 ? -half : half) * axis[static_cast<Eigen::Index>(c)];
  Method method;
  const BasisSet basis_a = BasisSet::build(aligned, BasisLibrary::load("6-31G*"));
  const BasisSet basis_t = BasisSet::build(tilted, BasisLibrary::load("6-31G*"));
  DirectJK jk_a(basis_a), jk_t(basis_t);
  const ScfResult scf_a = run_scf(aligned, basis_a, method, jk_a, tight()), scf_t = run_scf(tilted, basis_t, method, jk_t, tight());
  const RealTimePropagator prop_a(aligned, basis_a, method, scf_a, jk_a), prop_t(tilted, basis_t, method, scf_t, jk_t);

  // optical Kerr effect
  HarmonicFieldOptions harmonic;
  harmonic.omega = 0.1;
  harmonic.amplitude = 0.01;
  harmonic.ramp_time = 300.0;
  harmonic.n_periods = 2;
  harmonic.propagation.time_step = 0.2;
  harmonic.propagation.fock_tolerance = 1e-12;
  harmonic.optical = {0.0, 0.0, 1.0};
  const double gamma_bond = 8.0 * harmonic_response(prop_a, harmonic).cosine[2][1][2];
  const KerrAverage kerr = kerr_average(prop_t, harmonic);
  std::printf("non-linear, H2 6-31G*: gamma(-w; w, w, -w) along the bond %.4f, one fifth %.4f, average of the tilted molecule %.4f\n",
              gamma_bond, gamma_bond / 5.0, kerr.gamma);
  REQUIRE(kerr.steps_converged);
  CHECK_THAT(kerr.gamma, WithinRel(gamma_bond / 5.0, 1e-2));

  // Two-photon strength of the dark state of H2 in cc-pVDZ (half its energy is below the first allowed
  // excitation, as the method requires): the average must be the same for both orientations.
  const BasisSet big_a = BasisSet::build(aligned, BasisLibrary::load("cc-pVDZ"));
  const BasisSet big_t = BasisSet::build(tilted, BasisLibrary::load("cc-pVDZ"));
  DirectJK big_jk_a(big_a), big_jk_t(big_t);
  const ScfResult big_scf_a = run_scf(aligned, big_a, method, big_jk_a, tight()),
                  big_scf_t = run_scf(tilted, big_t, method, big_jk_t, tight());
  const RealTimePropagator big_prop_a(aligned, big_a, method, big_scf_a, big_jk_a),
      big_prop_t(tilted, big_t, method, big_scf_t, big_jk_t);
  const ResponseOperator op(aligned, big_a, method, big_scf_a, big_jk_a);
  ExcitationOptions excitation_options;
  excitation_options.n_states = 3;
  const ExcitationResult states = excitations(op, excitation_options);
  REQUIRE(states.converged);
  double dark = 0.0, first_allowed = 0.0;
  for (const Excitation& s : states.states) {
    if (dark == 0.0 && s.oscillator_strength < 1e-8) dark = s.energy;
    if (first_allowed == 0.0 && s.oscillator_strength > 1e-3) first_allowed = s.energy;
  }
  REQUIRE(dark > 0.0);
  REQUIRE(0.5 * dark < first_allowed - 0.1);
  TwoPhotonStrengthOptions options;
  options.state_energies = {dark};
  options.duration = 200.0;
  options.propagation.propagator = Propagator::Composition4;
  options.propagation.time_step = 0.2;
  options.propagation.fock_tolerance = 1e-11;
  options.propagation.max_fock_builds_per_step = 20;
  const TwoPhotonStrengths along_bond = two_photon_strengths(big_prop_a, options);
  const TwoPhotonStrengths averaged = two_photon_strengths(big_prop_t, options);
  REQUIRE(averaged.steps_converged);
  const double peak = two_photon_cross_section_gm(averaged.states, 0.5 * dark, 0.1 / 27.211386245988);
  std::printf("non-linear, H2 cc-pVDZ: dark state %.4f, averaged two-photon strength %.4f (aligned) and %.4f (tilted); "
              "cross-section at the peak with a half width of 0.1 eV: %.3f GM; %d Fock builds\n",
              dark, along_bond.states[0].strength, averaged.states[0].strength, peak, along_bond.fock_builds + averaged.fock_builds);
  CHECK(along_bond.states[0].strength > 0.0);
  CHECK_THAT(averaged.states[0].strength, WithinRel(along_bond.states[0].strength, 1e-2));
  // the line shape integrates to one: the peak is strength / (pi * half width) times the prefactor
  const double w = 0.5 * dark, hw = 0.1 / 27.211386245988;
  const double prefactor = 4.0 * std::pow(M_PI, 3) * w * w / (137.036 * 137.036) * std::pow(0.529177210903e-8, 4) * 2.4188843265864e-17 / 1e-50;
  CHECK_THAT(peak, WithinRel(prefactor * averaged.states[0].strength / (M_PI * hw), 1e-12));
}
