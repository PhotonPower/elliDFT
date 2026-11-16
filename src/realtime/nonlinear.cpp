#include "realtime/nonlinear.hpp"

#include <cmath>
#include <stdexcept>

#include "core/units.hpp"
#include "realtime/spectrum.hpp"

namespace ellidft {

namespace {

// Atomic unit of time in seconds (Wikipedia, "Hartree atomic units").
constexpr double atomic_time_s = 2.4188843265864e-17;
// Speed of light in atomic units, 1 / alpha. The page above gives "c = 1/alpha
// a.u. ~ 137 a.u."; the further digits are the commonly quoted value and were
// not taken from a consulted source.
constexpr double speed_of_light = 137.036;

}  // namespace

double smooth_switch(double x) {
  if (x <= 0.0) return 0.0;
  if (x >= 1.0) return 1.0;
  return x * x * x * x * (35.0 - 84.0 * x + 70.0 * x * x - 20.0 * x * x * x);
}

HarmonicResponse harmonic_response(const RealTimePropagator& propagator, const HarmonicFieldOptions& options) {
  if (!(options.omega > 0.0)) throw std::invalid_argument("harmonic_response: omega must be positive");
  if (!(options.ramp_time > 0.0) || options.n_periods < 1)
    throw std::invalid_argument("harmonic_response: needs a ramp time and at least one period");
  const double w = options.omega, ramp = options.ramp_time;
  const double period = 2.0 * M_PI / w;
  RealTimeOptions propagation = options.propagation;
  const double dt = propagation.time_step;
  const auto first = static_cast<std::size_t>(std::ceil(ramp / dt));
  propagation.n_steps = static_cast<int>(first) + static_cast<int>(std::ceil(options.n_periods * period / dt));

  const Vec3 f1 = options.optical, f0 = options.static_field;
  const DipoleOrders orders = dipole_orders(
      propagator,
      [=](double t) {
        const double s = smooth_switch(t / ramp), c = std::cos(w * t);
        return Vec3{s * (f0[0] + f1[0] * c), s * (f0[1] + f1[1] * c), s * (f0[2] + f1[2] * c)};
      },
      options.amplitude, propagation);

  HarmonicResponse result;
  result.fock_builds = orders.fock_builds;
  result.steps_converged = orders.steps_converged;
  result.total_time = dt * propagation.n_steps;

  // Least squares in the window: columns 1, cos(m w t), sin(m w t), m = 1..3.
  const std::size_t n_points = orders.first.size() - first;
  Eigen::MatrixXd basis(static_cast<Eigen::Index>(n_points), 7);
  for (std::size_t k = 0; k < n_points; ++k) {
    const double t = dt * static_cast<double>(first + k);
    const auto row = static_cast<Eigen::Index>(k);
    basis(row, 0) = 1.0;
    for (int m = 1; m <= 3; ++m) {
      basis(row, m) = std::cos(m * w * t);
      basis(row, 3 + m) = std::sin(m * w * t);
    }
  }
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(basis);
  const std::array<const std::vector<Vec3>*, 3> signals = {&orders.first, &orders.second, &orders.third};
  for (std::size_t n = 0; n < 3; ++n) {
    double signal_square = 0.0, residual_square = 0.0;
    for (std::size_t component = 0; component < 3; ++component) {
      Eigen::VectorXd y(static_cast<Eigen::Index>(n_points));
      for (std::size_t k = 0; k < n_points; ++k) y[static_cast<Eigen::Index>(k)] = (*signals[n])[first + k][component];
      const Eigen::VectorXd c = qr.solve(y);
      result.cosine[n][0][component] = c[0];
      for (int m = 1; m <= 3; ++m) {
        result.cosine[n][static_cast<std::size_t>(m)][component] = c[m];
        result.sine[n][static_cast<std::size_t>(m)][component] = c[3 + m];
      }
      signal_square += y.squaredNorm();
      residual_square += (y - basis * c).squaredNorm();
    }
    result.residual[n] = signal_square > 0.0 ? std::sqrt(residual_square / signal_square) : 0.0;
    result.signal_rms[n] = std::sqrt(signal_square / static_cast<double>(3 * n_points));
  }
  return result;
}

TwoPhotonResult two_photon_absorption(const RealTimePropagator& propagator, const TwoPhotonOptions& options) {
  if (!(options.omega > 0.0) || !(options.duration > 0.0) || !(options.amplitude > 0.0))
    throw std::invalid_argument("two_photon_absorption: omega, duration and amplitude must be positive");
  const double w = options.omega, T = options.duration, a = options.amplitude;
  double norm = 0.0;
  for (double x : options.polarization) norm += x * x;
  if (!(norm > 0.0)) throw std::invalid_argument("two_photon_absorption: no polarisation direction");
  const double scale = 1.0 / std::sqrt(norm);
  const Vec3 e{options.polarization[0] * scale, options.polarization[1] * scale, options.polarization[2] * scale};

  TwoPhotonResult result;
  auto absorbed = [&](double F0) {
    RealTimeOptions propagation = options.propagation;
    propagation.n_steps = static_cast<int>(std::ceil(T / propagation.time_step));
    const Vec3 constant = options.static_field;
    propagation.field = [=](double t) {
      if (t <= 0.0 || t >= T) return constant;
      const double s = std::sin(M_PI * t / T);
      const double f = F0 * s * s * std::cos(w * (t - 0.5 * T));
      return Vec3{constant[0] + f * e[0], constant[1] + f * e[1], constant[2] + f * e[2]};
    };
    const RealTimeTrajectory t = propagator.propagate({0.0, 0.0, 0.0}, propagation);
    result.fock_builds += t.fock_builds;
    result.steps_converged = result.steps_converged && t.steps_converged;
    // The trajectory holds the energy without the interaction with the
    // field. With a static field the conserved energy after the pulse is
    // E - mu.F_static.
    double change = t.energy.back() - t.energy.front();
    for (std::size_t k = 0; k < 3; ++k) change -= constant[k] * (t.dipole.back()[k] - t.dipole.front()[k]);
    return change;
  };
  const double W1 = absorbed(a), W2 = absorbed(2.0 * a);
  result.one_photon_coefficient = (16.0 * W1 - W2) / (12.0 * a * a);
  result.two_photon_coefficient = (W2 - 4.0 * W1) / (12.0 * a * a * a * a);

  const double envelope = T * 35.0 / 128.0;   // integral of sin^8(pi t / T) over the pulse
  const double factor = 8.0 * M_PI / speed_of_light;
  result.cross_section_au = result.two_photon_coefficient * w * factor * factor / envelope;
  // Bohr^4 * (atomic unit of time) in cm^4 s, in units of 1e-50 (1 GM =
  // 1e-50 cm^4 s / photon: Wikipedia, "Two-photon absorption").
  const double bohr_cm = units::bohr_in_angstrom * 1e-8;
  result.cross_section_gm = result.cross_section_au * bohr_cm * bohr_cm * bohr_cm * bohr_cm * atomic_time_s / 1e-50;
  return result;
}

std::array<AveragingDirection, 7> averaging_directions() {
  const double d = 1.0 / std::sqrt(3.0), axis = 2.0 / 15.0, diagonal = 3.0 / 20.0;
  return {{{{1.0, 0.0, 0.0}, axis}, {{0.0, 1.0, 0.0}, axis}, {{0.0, 0.0, 1.0}, axis},
           {{d, d, d}, diagonal}, {{d, d, -d}, diagonal}, {{d, -d, d}, diagonal}, {{-d, d, d}, diagonal}}};
}

KerrAverage kerr_average(const RealTimePropagator& propagator, const HarmonicFieldOptions& options) {
  KerrAverage result;
  const auto directions = averaging_directions();
  double largest_left = 0.0, largest_signal = 0.0;
  for (std::size_t k = 0; k < directions.size(); ++k) {
    HarmonicFieldOptions one = options;
    one.optical = directions[k].u;
    one.static_field = {0.0, 0.0, 0.0};
    const HarmonicResponse r = harmonic_response(propagator, one);
    const Vec3& u = directions[k].u;
    // third order, first harmonic: 1/8 gamma(-w; w, w, -w)
    const double in_phase = 8.0 * (r.cosine[2][1][0] * u[0] + r.cosine[2][1][1] * u[1] + r.cosine[2][1][2] * u[2]);
    const double quadrature = 8.0 * (r.sine[2][1][0] * u[0] + r.sine[2][1][1] * u[1] + r.sine[2][1][2] * u[2]);
    result.by_direction[k] = in_phase;
    result.gamma += directions[k].weight * in_phase;
    result.gamma_quadrature += directions[k].weight * quadrature;
    largest_left = std::max(largest_left, r.residual[2] * r.signal_rms[2]);
    largest_signal = std::max(largest_signal, r.signal_rms[2]);
    result.fock_builds += r.fock_builds;
    result.steps_converged = result.steps_converged && r.steps_converged;
  }
  result.largest_residual = largest_signal > 0.0 ? largest_left / largest_signal : 0.0;
  return result;
}

NonlinearIndex nonlinear_index(double gamma_au, double number_density_per_A3, double n0) {
  if (!(number_density_per_A3 > 0.0) || !(n0 >= 1.0))
    throw std::invalid_argument("nonlinear_index: needs a positive number density and n0 >= 1");
  NonlinearIndex result;
  const double bohr = units::bohr_in_angstrom;
  const double N = number_density_per_A3 * bohr * bohr * bohr;   // per Bohr^3
  const double L = (n0 * n0 + 2.0) / 3.0;
  result.local_field_factor = L;
  result.n2_au = 2.0 * M_PI * M_PI * N * L * L * L * L * gamma_au / (n0 * n0 * speed_of_light);
  // Atomic unit of intensity: Hartree / (atomic unit of time * Bohr^2), with
  // the Hartree in Joule from its value in eV and the elementary charge
  // 1.602176634e-19 C (exact in the SI).
  const double hartree_J = units::hartree_in_ev * 1.602176634e-19;
  const double bohr_m = bohr * 1e-10;
  const double intensity_W_per_m2 = hartree_J / (atomic_time_s * bohr_m * bohr_m);
  result.n2_m2_per_W = result.n2_au / intensity_W_per_m2;
  result.n2_cm2_per_W = result.n2_m2_per_W * 1e4;
  return result;
}

namespace {

// Simpson's rule for a complex integrand on [0, T]; the integrands here are
// smooth and vanish at both ends.
template <typename Function>
std::complex<double> integrate(const Function& f, double T, double fastest_frequency) {
  // at least 64 points per period of the fastest oscillation
  int n = static_cast<int>(std::ceil(64.0 * T * std::max(fastest_frequency, 2.0 * M_PI / T) / (2.0 * M_PI)));
  n = std::max(n + n % 2, 256);
  const double h = T / n;
  std::complex<double> sum = f(0.0) + f(T);
  for (int k = 1; k < n; ++k) sum += (k % 2 == 1 ? 4.0 : 2.0) * f(h * k);
  return sum * h / 3.0;
}

}  // namespace

std::complex<double> pulse_transform(double omega0, double duration, double omega) {
  const double T = duration;
  return integrate(
      [=](double t) {
        const double s = std::sin(M_PI * t / T);
        return s * s * std::cos(omega0 * (t - 0.5 * T)) * std::exp(std::complex<double>(0.0, omega * t));
      },
      T, std::abs(omega0) + std::abs(omega));
}

std::complex<double> envelope_squared_transform(double duration, double delta) {
  const double T = duration;
  return integrate(
      [=](double t) {
        const double s = std::sin(M_PI * t / T);
        return s * s * s * s * std::exp(std::complex<double>(0.0, delta * t));
      },
      T, std::abs(delta));
}

TwoPhotonStrengths two_photon_strengths(const RealTimePropagator& propagator, const TwoPhotonStrengthOptions& options) {
  const std::size_t n_states = options.state_energies.size();
  if (n_states == 0) throw std::invalid_argument("two_photon_strengths: no states");
  for (double energy : options.state_energies)
    if (!(energy > 0.0)) throw std::invalid_argument("two_photon_strengths: excitation energies must be positive");

  // W4_k = sum_f M_kf strength_f for pulses at w_k = w_k(state) / 2
  Eigen::MatrixXd M(static_cast<Eigen::Index>(n_states), static_cast<Eigen::Index>(n_states));
  for (std::size_t k = 0; k < n_states; ++k)
    for (std::size_t f = 0; f < n_states; ++f) {
      const double delta = options.state_energies[f] - options.state_energies[k];
      M(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(f)) =
          options.state_energies[f] * std::norm(envelope_squared_transform(options.duration, delta)) / 64.0;
    }
  const Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver(M);

  TwoPhotonStrengths result;
  result.states.resize(n_states);
  for (std::size_t f = 0; f < n_states; ++f) result.states[f].energy = options.state_energies[f];

  auto add_direction = [&](const Vec3& u, double weight) {
    Eigen::VectorXd W4(static_cast<Eigen::Index>(n_states));
    for (std::size_t k = 0; k < n_states; ++k) {
      TwoPhotonOptions pulse;
      pulse.omega = 0.5 * options.state_energies[k];
      pulse.polarization = u;
      pulse.amplitude = options.amplitude;
      pulse.duration = options.duration;
      pulse.propagation = options.propagation;
      const TwoPhotonResult r = two_photon_absorption(propagator, pulse);
      W4[static_cast<Eigen::Index>(k)] = r.two_photon_coefficient;
      result.fock_builds += r.fock_builds;
      result.steps_converged = result.steps_converged && r.steps_converged;
    }
    const Eigen::VectorXd strengths = solver.solve(W4);
    for (std::size_t f = 0; f < n_states; ++f) result.states[f].strength += weight * strengths[static_cast<Eigen::Index>(f)];
  };
  if (options.average) {
    for (const AveragingDirection& d : averaging_directions()) add_direction(d.u, d.weight);
  } else {
    add_direction(options.polarization, 1.0);
  }
  return result;
}

double two_photon_cross_section_gm(const std::vector<TwoPhotonState>& states, double omega, double half_width) {
  if (!(half_width > 0.0)) throw std::invalid_argument("two_photon_cross_section_gm: the line width must be positive");
  double sum = 0.0;
  for (const TwoPhotonState& s : states) {
    const double x = 2.0 * omega - s.energy;
    sum += s.strength * half_width / (M_PI * (x * x + half_width * half_width));
  }
  const double alpha = 1.0 / speed_of_light;
  const double sigma_au = 4.0 * M_PI * M_PI * M_PI * alpha * alpha * omega * omega * sum;
  const double bohr_cm = units::bohr_in_angstrom * 1e-8;
  return sigma_au * bohr_cm * bohr_cm * bohr_cm * bohr_cm * atomic_time_s / 1e-50;
}

}  // namespace ellidft
