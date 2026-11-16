// Module realtime: propagation of the orbitals after a field pulse.
// No reference program. The central check is against the module response:
// for a weak pulse the damped Fourier transform of the dipole moment must be
// the polarisability alpha(omega + i gamma) of linear response theory, which
// is computed there in a completely different way (frequency domain,
// exchange-correlation kernel instead of the potential).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <complex>
#include <cstdio>
#include <string>

#include "core/xyz.hpp"
#include "integrals/one_electron.hpp"
#include "realtime/propagator.hpp"
#include "realtime/spectrum.hpp"
#include "response/polarizability.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Molecule h2() { return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2.xyz"); }

// Water in a general orientation: no tensor component vanishes by symmetry.
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

// Small grids: the real-time and the linear response route use the same
// one, so their agreement does not depend on its quality.
Method method_of(const std::string& functional, const GridSpec& grid = {30, 110}) {
  Method method;
  method.functional = Functional::by_name(functional);
  method.grid = grid;
  return method;
}

// Largest deviation of the real-time polarisability from linear response,
// relative to the largest element of the tensor, over some frequencies.
struct Comparison {
  double deviation = 0.0;
  DipoleResponse response;
};
Comparison compare_with_linear_response(const Molecule& molecule, const std::string& basis_name,
                                        const std::string& functional, const RealTimeOptions& options,
                                        const KickOptions& kick, double gamma, const GridSpec& grid = {30, 110}) {
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(basis_name));
  const Method method = method_of(functional, grid);
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  REQUIRE(scf.converged);

  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  Comparison result;
  result.response = dipole_response(propagator, options, kick);

  const ResponseOperator op(molecule, basis, method, scf, jk);
  LinearResponseOptions tight_response;
  tight_response.tolerance = 1e-11;
  for (double omega : {0.0, 0.08, 0.2, 0.45}) {
    const Polarizability reference = polarizability(op, omega, gamma, tight_response);
    REQUIRE(reference.converged);
    const Eigen::Matrix3cd alpha = polarizability(result.response, omega, gamma);
    double scale = 0.0, difference = 0.0;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        if (!kick.directions[static_cast<std::size_t>(j)]) continue;
        const std::complex<double> ref(reference.real(i, j), reference.imag(i, j));
        scale = std::max(scale, std::abs(ref));
        difference = std::max(difference, std::abs(alpha(i, j) - ref));
      }
    result.deviation = std::max(result.deviation, difference / scale);
  }
  std::printf("real time vs linear response: %-10s %-8s %-9s dt = %.3f  T = %.0f  relative deviation %.2e  "
              "(%d Fock builds, energy drift %.1e)\n",
              molecule.name.c_str(), basis_name.c_str(), functional.c_str(), options.time_step,
              result.response.total_time(), result.deviation, result.response.fock_builds, result.response.energy_drift);
  std::fflush(stdout);
  return result;
}

}  // namespace

TEST_CASE("real time: Fock matrix and energy of the ground state are those of the SCF", "[realtime]") {
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("cc-pVDZ"));
  for (const char* functional : {"HF", "PBE", "PBE0", "CAM-B3LYP"}) {
    INFO(functional);
    const Method method = method_of(functional);
    DirectJK jk(basis);
    const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
    REQUIRE(scf.converged);
    const RealTimePropagator propagator(molecule, basis, method, scf, jk);
    const auto result = propagator.fock(scf.density.cast<std::complex<double>>());
    CHECK_THAT(result.energy, WithinAbs(scf.energy_total, 1e-10));
    CHECK((result.fock.real() - scf.fock).cwiseAbs().maxCoeff() < 1e-9);
    CHECK(result.fock.imag().cwiseAbs().maxCoeff() == 0.0);
    CHECK(propagator.n_occupied() == 1);
    CHECK(propagator.n_orbitals() == scf.n_orbitals);
  }
}

TEST_CASE("real time: without a pulse nothing moves", "[realtime]") {
  const Molecule molecule = tilted_water();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  const Method method = method_of("HF");
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  RealTimeOptions options;
  options.n_steps = 50;
  const RealTimeTrajectory t = propagator.propagate({0.0, 0.0, 0.0}, options);
  REQUIRE(t.dipole.size() == 51);
  for (std::size_t n = 0; n < t.dipole.size(); ++n) {
    for (std::size_t k = 0; k < 3; ++k) CHECK_THAT(t.dipole[n][k], WithinAbs(scf.dipole[k], 1e-9));
    CHECK_THAT(t.energy[n], WithinAbs(scf.energy_total, 1e-10));
  }
  CHECK(t.steps_converged);
  CHECK(t.max_fock_builds_in_a_step == 1);   // the extrapolated midpoint is already right
}

TEST_CASE("real time: the Fock matrix of a complex density is Hermitian", "[realtime]") {
  const Molecule molecule = tilted_water();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  for (const char* functional : {"HF", "PBE", "PBE0", "CAM-B3LYP"}) {
    INFO(functional);
    const Method method = method_of(functional);
    DirectJK jk(basis);
    const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
    const RealTimePropagator propagator(molecule, basis, method, scf, jk);

    // a complex density: ground-state orbitals with an admixture i * 0.2 * x
    const Eigen::MatrixXd D = dipole_matrices(basis)[0];
    Eigen::MatrixXcd c = scf.coefficients.leftCols(scf.n_occupied).cast<std::complex<double>>();
    c = c + std::complex<double>(0.0, 0.2) * (scf.overlap.inverse() * D * c);
    const auto fe = propagator.fock(2.0 * c * c.adjoint());
    CHECK((fe.fock - fe.fock.adjoint()).cwiseAbs().maxCoeff() < 1e-12);
    // only exact exchange sees the imaginary part of the density
    if (method.functional.exact_exchange != 0.0 || method.functional.long_range_exchange != 0.0)
      CHECK(fe.fock.imag().cwiseAbs().maxCoeff() > 1e-3);
    else
      CHECK(fe.fock.imag().cwiseAbs().maxCoeff() == 0.0);
  }
}

TEST_CASE("real time: energy and orthonormality after a strong pulse", "[realtime]") {
  // After the pulse the Hamiltonian no longer depends on time: the energy is
  // conserved by the exact dynamics. The midpoint propagator and the
  // composition built from it conserve it to rounding at any step (Hartree-
  // Fock). The Magnus propagator of this module does not (it
  // uses Fock matrices of earlier steps); its error must fall at least like
  // the fourth power of the step once the step resolves the fastest motion.
  const Molecule molecule = tilted_water();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  const Method method = method_of("HF");
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const Vec3 pulse{0.03, -0.02, 0.05};

  auto run = [&](Propagator scheme, double dt) {
    RealTimeOptions options;
    options.propagator = scheme;
    options.time_step = dt;
    options.n_steps = static_cast<int>(10.0 / dt + 0.5);
    options.fock_tolerance = 1e-11;
    options.max_fock_builds_per_step = 30;
    const RealTimeTrajectory t = propagator.propagate(pulse, options);
    CHECK(t.steps_converged);
    CHECK(t.orthonormality_error < 1e-11);
    CHECK(t.energy.front() - scf.energy_total > 1e-3);   // the pulse puts energy in
    return t.energy_drift;
  };
  const double midpoint = run(Propagator::Midpoint, 0.1), composition = run(Propagator::Composition4, 0.1);
  const double magnus_coarse = run(Propagator::Magnus4, 0.025), magnus_fine = run(Propagator::Magnus4, 0.0125);
  std::printf("real time, strong pulse: energy drift midpoint %.1e and Composition4 %.1e (dt 0.1), Magnus4 %.1e (dt 0.025) and %.1e (dt 0.0125)\n",
              midpoint, composition, magnus_coarse, magnus_fine);
  CHECK(midpoint < 1e-10);
  CHECK(composition < 1e-10);
  CHECK(magnus_fine < 1e-6);
  CHECK(magnus_coarse / magnus_fine > 16.0);
}

TEST_CASE("real time: transform of a known signal", "[realtime]") {
  // d mu / d kappa = c sin(W t)  ->  alpha(z) = c W / (W^2 - z^2), z = omega + i gamma
  // (integral_0^inf sin(W t) e^{i z t} dt). Two lines with different weights.
  const double W1 = 0.41, c1 = 3.0, W2 = 1.7, c2 = 0.8, gamma = 0.05;
  DipoleResponse response;
  response.time_step = 0.05;
  for (int n = 0; n <= 6001; ++n) {   // odd number of intervals: last three by the 3/8 rule
    const double t = response.time_step * n;
    Eigen::Matrix3d r = Eigen::Matrix3d::Zero();
    r(0, 0) = c1 * std::sin(W1 * t) + c2 * std::sin(W2 * t);
    r(1, 2) = c2 * std::sin(W2 * t);
    response.response.push_back(r);
  }
  for (double omega : {0.0, 0.2, 0.41, 1.0, 1.7}) {
    const std::complex<double> z(omega, gamma);
    const std::complex<double> a1 = c1 * W1 / (W1 * W1 - z * z), a2 = c2 * W2 / (W2 * W2 - z * z);
    const Eigen::Matrix3cd alpha = polarizability(response, omega, gamma);
    // Simpson's rule: error of the order (dt * frequency)^4, here below 1e-5
    CHECK_THAT(std::abs(alpha(0, 0) - (a1 + a2)), WithinAbs(0.0, 1e-5 * std::abs(a1 + a2)));
    CHECK_THAT(std::abs(alpha(1, 2) - a2), WithinAbs(0.0, 1e-5 * std::abs(a2)));
    CHECK(std::abs(alpha(2, 1)) == 0.0);
    CHECK(alpha(0, 0).imag() >= 0.0);   // absorption
  }
  CHECK_THROWS_AS(polarizability(response, 0.1, 0.0), std::invalid_argument);
  response.response.resize(2);
  CHECK_THROWS_AS(polarizability(response, 0.1, 0.1), std::invalid_argument);
}

// Criterion proposed to the client: relative deviation from linear response
// at most 1e-4 for H2 and H2O with HF, PBE and PBE0; electron number
// (orthonormality of the orbitals) within 1e-10; energy after the pulse
// within 1e-7 Hartree.
//
// Short versions for every test run: H2, strong damping (short signal), for
// the density functionals a smaller grid, one pulse direction and one
// propagation (H2 has a centre of symmetry: no even orders in the field).
TEST_CASE("real time: polarisability of H2 equals linear response", "[realtime]") {
  RealTimeOptions options;
  options.time_step = 0.2;
  options.n_steps = 225;        // T = 45, exp(-gamma T) = 1.4e-6
  const double gamma = 0.3;
  SECTION("Hartree-Fock, perpendicular and parallel to the bond") {
    KickOptions kick;
    kick.directions = {true, false, true};
    const Comparison c = compare_with_linear_response(h2(), "cc-pVDZ", "HF", options, kick, gamma);
    CHECK(c.deviation <= 1e-4);
    CHECK(c.response.steps_converged);
    CHECK(c.response.propagations == 4);
    CHECK(c.response.orthonormality_error < 1e-10);
    CHECK(c.response.energy_drift < 1e-7);
  }
  SECTION("PBE and PBE0, along the bond") {
    KickOptions kick;
    kick.directions = {false, false, true};
    kick.two_sided = false;
    for (const char* functional : {"PBE", "PBE0"}) {
      INFO(functional);
      const Comparison c = compare_with_linear_response(h2(), "cc-pVDZ", functional, options, kick, gamma, {20, 50});
      CHECK(c.deviation <= 1e-4);
      CHECK(c.response.steps_converged);
      CHECK(c.response.orthonormality_error < 1e-10);
      CHECK(c.response.energy_drift < 1e-7);
    }
  }
}

TEST_CASE("real time: polarisability of H2 equals linear response (all functionals, both directions)",
          "[realtime][slow]") {
  RealTimeOptions options;
  options.time_step = 0.2;
  options.n_steps = 450;        // T = 90, exp(-gamma T) = 1.4e-6
  KickOptions kick;
  kick.directions = {true, false, true};
  for (const char* functional : {"HF", "PBE", "PBE0", "CAM-B3LYP"}) {
    INFO(functional);
    const Comparison c = compare_with_linear_response(h2(), "cc-pVDZ", functional, options, kick, 0.15);
    CHECK(c.deviation <= 1e-4);
    CHECK(c.response.steps_converged);
    CHECK(c.response.orthonormality_error < 1e-10);
    CHECK(c.response.energy_drift < 1e-7);
  }
}

TEST_CASE("real time: full polarisability tensor of water equals linear response", "[realtime]") {
  // Water has no centre of symmetry: pulses of both signs. The 1s orbital of
  // oxygen brings frequencies of 20 Hartree; hence the shorter time step.
  RealTimeOptions options;
  options.time_step = 0.1;
  options.n_steps = 450;        // T = 45, exp(-gamma T) = 1.4e-6
  const Comparison c = compare_with_linear_response(tilted_water(), "STO-3G", "HF", options, KickOptions{}, 0.3);
  CHECK(c.deviation <= 1e-4);
  CHECK(c.response.propagations == 6);
  CHECK(c.response.steps_converged);
  CHECK(c.response.orthonormality_error < 1e-10);
  CHECK(c.response.energy_drift < 1e-7);
}

TEST_CASE("real time: water with PBE and PBE0 against linear response", "[realtime][slow]") {
  RealTimeOptions options;
  options.time_step = 0.1;
  options.n_steps = 900;
  for (const char* functional : {"PBE", "PBE0"}) {
    INFO(functional);
    const Comparison c = compare_with_linear_response(tilted_water(), "STO-3G", functional, options, KickOptions{}, 0.15);
    CHECK(c.deviation <= 1e-4);
    CHECK(c.response.steps_converged);
    CHECK(c.response.energy_drift < 1e-7);
  }
}

TEST_CASE("real time: orders of convergence of the propagators", "[realtime]") {
  // Halving the time step must reduce the error by about 4 (midpoint) and
  // by about 16 (Magnus and composition, fourth order). Damping 0.4 and
  // T = 40: the neglected tail, exp(-16) = 1e-7, is below all errors
  // compared here.
  KickOptions kick;
  kick.directions = {false, false, true};
  auto ratio = [&](Propagator scheme) {
    double error[2];
    int k = 0;
    for (double dt : {0.4, 0.2}) {
      RealTimeOptions options;
      options.propagator = scheme;
      options.time_step = dt;
      options.n_steps = static_cast<int>(40.0 / dt + 0.5);
      options.fock_tolerance = 1e-11;
      error[k++] = compare_with_linear_response(h2(), "cc-pVDZ", "HF", options, kick, 0.4).deviation;
    }
    return error[0] / error[1];
  };
  const double midpoint = ratio(Propagator::Midpoint), magnus = ratio(Propagator::Magnus4),
               composition = ratio(Propagator::Composition4);
  std::printf("real time: error ratio for half the time step: midpoint %.2f, Magnus4 %.2f, Composition4 %.2f\n",
              midpoint, magnus, composition);
  CHECK(midpoint > 3.5);
  CHECK(midpoint < 4.5);
  CHECK(magnus > 10.0);
  CHECK(composition > 10.0);
}

namespace {

// Smooth switch from 0 to 1 on [0, 1] whose first three derivatives vanish
// at both ends (polynomial of degree 7).
double smooth_step(double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x * x * x * x * (35.0 - 84.0 * x + 70.0 * x * x - 20.0 * x * x * x);
}

}  // namespace

TEST_CASE("real time: the energy absorbed from a pulse equals the work done by the field", "[realtime]") {
  // For a field along one axis, d E0 / dt = F(t) d mu / dt, with E0 the
  // energy without the interaction (from d<H0>/dt = i <[H, H0]>). For a
  // pulse that starts and ends at zero: E0(T) - E0(0) = - integral F'(t) mu(t) dt.
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("cc-pVDZ"));
  const Method method = method_of("HF");
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const double F0 = 0.02, w = 0.3, Tp = 40.0;
  auto field = [=](double t) {
    return (t < 0.0 || t > Tp) ? 0.0 : F0 * std::pow(std::sin(M_PI * t / Tp), 2) * std::cos(w * t);
  };
  auto field_derivative = [=](double t) {
    if (t < 0.0 || t > Tp) return 0.0;
    const double s = std::sin(M_PI * t / Tp), c = std::cos(M_PI * t / Tp);
    return F0 * (2.0 * s * c * M_PI / Tp * std::cos(w * t) - s * s * w * std::sin(w * t));
  };
  struct Case {
    Propagator scheme;
    double dt, tolerance;
  };
  // tolerances: one order of magnitude above the measured differences
  for (const Case& c : {Case{Propagator::Midpoint, 0.1, 3e-6}, Case{Propagator::Magnus4, 0.1, 2e-8},
                        Case{Propagator::Composition4, 0.2, 3e-7}}) {
    RealTimeOptions options;
    options.propagator = c.scheme;
    options.time_step = c.dt;
    options.n_steps = static_cast<int>(Tp / c.dt + 0.5);
    options.field = [&](double t) { return Vec3{0.0, 0.0, field(t)}; };
    const RealTimeTrajectory tr = propagator.propagate({0.0, 0.0, 0.0}, options);
    REQUIRE(tr.steps_converged);
    double work = 0.0;   // Simpson's rule
    const std::size_t N = tr.dipole.size() - 1;
    REQUIRE(N % 2 == 0);
    for (std::size_t k = 0; k <= N; ++k) {
      const double weight = (k == 0 || k == N) ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0);
      work -= weight * c.dt / 3.0 * field_derivative(c.dt * static_cast<double>(k)) * tr.dipole[k][2];
    }
    const double absorbed = tr.energy.back() - tr.energy.front();
    std::printf("real time, pulse: propagator %d, dt %.1f: absorbed %.6e Hartree, work %.6e, difference %.1e\n",
                static_cast<int>(c.scheme), c.dt, absorbed, work, absorbed - work);
    CHECK(absorbed > 1e-3);
    CHECK_THAT(absorbed, WithinAbs(work, c.tolerance));
  }
}

TEST_CASE("real time: a slowly switched field gives the dipole moment of the SCF in that field", "[realtime]") {
  // Strong static field: the response is visibly non-linear. The SCF with
  // the field in the Hamiltonian is an independent route to the same state.
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("cc-pVDZ"));
  const Method method = method_of("HF");
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const ResponseOperator op(molecule, basis, method, scf, jk);
  const double alpha = polarizability(op).real(2, 2);

  const double F = 0.1, ramp = 200.0, hold = 20.0, dt = 0.2;
  ScfOptions in_field = tight();
  in_field.electric_field = {0.0, 0.0, F};
  const ScfResult reference = run_scf(molecule, basis, method, jk, in_field);
  REQUIRE(reference.converged);
  const double linear = alpha * F;
  CHECK(std::abs(reference.dipole[2] - linear) > 5e-3 * linear);   // non-linear part above 0.5 %

  RealTimeOptions options;
  options.time_step = dt;
  options.n_steps = static_cast<int>((ramp + hold) / dt + 0.5);
  options.field = [=](double t) { return Vec3{0.0, 0.0, F * smooth_step(t / ramp)}; };
  const RealTimeTrajectory tr = propagator.propagate({0.0, 0.0, 0.0}, options);
  REQUIRE(tr.steps_converged);
  double deviation = 0.0;
  for (std::size_t k = static_cast<std::size_t>(ramp / dt + 0.5); k < tr.dipole.size(); ++k)
    deviation = std::max(deviation, std::abs(tr.dipole[k][2] - reference.dipole[2]));
  // E_SCF(F) = E0[P] - mu.F  ->  E0[P] = E_SCF(F) + mu.F
  const double energy_deviation = tr.energy.back() - (reference.energy_total + reference.dipole[2] * F);
  std::printf("real time, static field %.2f on H2: mu(SCF) = %.8f, alpha F = %.8f, largest deviation while held %.1e, energy %.1e\n",
              F, reference.dipole[2], linear, deviation, energy_deviation);
  CHECK(deviation < 1e-4 * reference.dipole[2]);
  CHECK(std::abs(energy_deviation) < 1e-6);
}

TEST_CASE("real time: dipole moment order by order in the field", "[realtime]") {
  // Water (no centre of symmetry), field switched on slowly along a general
  // direction u. In the static limit mu1 = alpha.u (linear response code);
  // mu2 and mu3 are compared with the same finite differences of SCF
  // calculations in static fields.
  const Molecule molecule = tilted_water();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  const Method method = method_of("HF");
  DirectJK jk(basis);
  const ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  const ResponseOperator op(molecule, basis, method, scf, jk);
  const Eigen::Matrix3d alpha = polarizability(op).real;
  const Eigen::Vector3d u = Eigen::Vector3d(1.0, -2.0, 1.5).normalized();
  const Eigen::Vector3d alpha_u = alpha * u;

  const double a = 0.02, ramp = 160.0, hold = 10.0, dt = 0.2;
  auto scf_dipole = [&](double lambda) {
    ScfOptions in_field = tight();
    in_field.electric_field = {lambda * u[0], lambda * u[1], lambda * u[2]};
    const ScfResult r = run_scf(molecule, basis, method, jk, in_field);
    REQUIRE(r.converged);
    return r.dipole;
  };
  const Vec3 p1 = scf_dipole(a), m1 = scf_dipole(-a), p2 = scf_dipole(2.0 * a), m2 = scf_dipole(-2.0 * a);

  RealTimeOptions options;
  options.time_step = dt;
  options.n_steps = static_cast<int>((ramp + hold) / dt + 0.5);
  options.fock_tolerance = 1e-11;   // the third order is a difference of dipole moments divided by 6 a^3
  const DipoleOrders orders = dipole_orders(
      propagator, [=](double t) { const double s = smooth_step(t / ramp); return Vec3{s * u[0], s * u[1], s * u[2]}; }, a,
      options);
  REQUIRE(orders.steps_converged);
  REQUIRE(orders.first.size() == static_cast<std::size_t>(options.n_steps) + 1);

  double scale1 = 0.0, scale2 = 0.0, scale3 = 0.0, dev1 = 0.0, dev2 = 0.0, dev3 = 0.0;
  for (std::size_t k = 0; k < 3; ++k) {
    const double o1 = 0.5 * (p1[k] - m1[k]), o2 = 0.5 * (p2[k] - m2[k]);
    const double e1 = 0.5 * (p1[k] + m1[k]) - scf.dipole[k], e2 = 0.5 * (p2[k] + m2[k]) - scf.dipole[k];
    const double mu2 = (16.0 * e1 - e2) / (12.0 * a * a), mu3 = (o2 - 2.0 * o1) / (6.0 * a * a * a);
    const Vec3& r1 = orders.first.back();
    const Vec3& r2 = orders.second.back();
    const Vec3& r3 = orders.third.back();
    scale1 = std::max(scale1, std::abs(alpha_u[static_cast<Eigen::Index>(k)]));
    scale2 = std::max(scale2, std::abs(mu2));
    scale3 = std::max(scale3, std::abs(mu3));
    dev1 = std::max(dev1, std::abs(r1[k] - alpha_u[static_cast<Eigen::Index>(k)]));
    dev2 = std::max(dev2, std::abs(r2[k] - mu2));
    dev3 = std::max(dev3, std::abs(r3[k] - mu3));
  }
  std::printf("real time, orders in the field (water): |mu1| %.4f dev %.1e, |mu2| %.4f dev %.1e, |mu3| %.3f dev %.1e (%d Fock builds)\n",
              scale1, dev1, scale2, dev2, scale3, dev3, orders.fock_builds);
  CHECK(scale2 > 0.1);   // the second order is there
  CHECK(dev1 < 1e-4 * scale1);
  CHECK(dev2 < 1e-3 * scale2);
  CHECK(dev3 < 1e-2 * scale3);
  // before the field is switched on nothing has happened
  for (std::size_t k = 0; k < 3; ++k) CHECK(orders.first.front()[k] == 0.0);
  CHECK_THROWS_AS(dipole_orders(propagator, nullptr, a, options), std::invalid_argument);
  CHECK_THROWS_AS(dipole_orders(propagator, [](double) { return Vec3{0.0, 0.0, 1.0}; }, 0.0, options), std::invalid_argument);
}

TEST_CASE("real time: invalid requests", "[realtime]") {
  const Molecule molecule = h2();
  const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load("STO-3G"));
  const Method method = method_of("HF");
  DirectJK jk(basis);
  ScfResult scf = run_scf(molecule, basis, method, jk, tight());
  const RealTimePropagator propagator(molecule, basis, method, scf, jk);
  RealTimeOptions options;
  options.time_step = 0.0;
  CHECK_THROWS_AS(propagator.propagate({0.0, 0.0, 1e-4}, options), std::invalid_argument);
  options = {};
  options.n_steps = 0;
  CHECK_THROWS_AS(propagator.propagate({0.0, 0.0, 1e-4}, options), std::invalid_argument);
  options = {};
  options.max_fock_builds_per_step = 0;
  CHECK_THROWS_AS(propagator.propagate({0.0, 0.0, 1e-4}, options), std::invalid_argument);
  KickOptions kick;
  kick.strength = 0.0;
  CHECK_THROWS_AS(dipole_response(propagator, RealTimeOptions{}, kick), std::invalid_argument);
  scf.converged = false;
  CHECK_THROWS_AS(RealTimePropagator(molecule, basis, method, scf, jk), std::invalid_argument);
}
