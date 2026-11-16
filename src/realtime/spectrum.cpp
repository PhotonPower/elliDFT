#include "realtime/spectrum.hpp"

#include <cmath>
#include <stdexcept>

namespace ellidft {

DipoleResponse dipole_response(const RealTimePropagator& propagator, const RealTimeOptions& options,
                               const KickOptions& kick) {
  if (!(kick.strength > 0.0)) throw std::invalid_argument("dipole_response: kick strength must be positive");
  DipoleResponse result;
  result.time_step = options.time_step;
  result.response.assign(static_cast<std::size_t>(options.n_steps) + 1, Eigen::Matrix3d::Zero());
  const Vec3 ground = propagator.ground_state_dipole();

  auto collect = [&](const RealTimeTrajectory& t) {
    ++result.propagations;
    result.fock_builds += t.fock_builds;
    result.steps_converged = result.steps_converged && t.steps_converged;
    result.energy_drift = std::max(result.energy_drift, t.energy_drift);
    result.orthonormality_error = std::max(result.orthonormality_error, t.orthonormality_error);
  };

  for (std::size_t j = 0; j < 3; ++j) {
    if (!kick.directions[j]) continue;
    Vec3 plus{0.0, 0.0, 0.0};
    plus[j] = kick.strength;
    const RealTimeTrajectory forward = propagator.propagate(plus, options);
    collect(forward);
    if (kick.two_sided) {
      Vec3 minus{0.0, 0.0, 0.0};
      minus[j] = -kick.strength;
      const RealTimeTrajectory backward = propagator.propagate(minus, options);
      collect(backward);
      for (std::size_t n = 0; n < result.response.size(); ++n)
        for (std::size_t i = 0; i < 3; ++i)
          result.response[n](static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
              (forward.dipole[n][i] - backward.dipole[n][i]) / (2.0 * kick.strength);
    } else {
      for (std::size_t n = 0; n < result.response.size(); ++n)
        for (std::size_t i = 0; i < 3; ++i)
          result.response[n](static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
              (forward.dipole[n][i] - ground[i]) / kick.strength;
    }
  }
  return result;
}

Eigen::Matrix3cd polarizability(const DipoleResponse& response, double omega, double gamma) {
  const std::size_t n_points = response.response.size();
  if (n_points < 3) throw std::invalid_argument("polarizability: needs at least two time steps");
  if (!(gamma > 0.0)) throw std::invalid_argument("polarizability: the damping gamma must be positive");
  if (!(response.time_step > 0.0)) throw std::invalid_argument("polarizability: time step must be positive");

  // Composite Simpson rule; with an odd number of intervals the last three
  // are integrated with the 3/8 rule.
  const std::size_t n_intervals = n_points - 1;
  std::vector<double> weights(n_points, 0.0);
  const std::size_t simpson_intervals = n_intervals % 2 == 0 ? n_intervals : n_intervals - 3;
  for (std::size_t k = 0; k + 2 <= simpson_intervals; k += 2) {
    weights[k] += 1.0 / 3.0;
    weights[k + 1] += 4.0 / 3.0;
    weights[k + 2] += 1.0 / 3.0;
  }
  if (simpson_intervals != n_intervals) {
    const std::size_t k = simpson_intervals;
    weights[k] += 3.0 / 8.0;
    weights[k + 1] += 9.0 / 8.0;
    weights[k + 2] += 9.0 / 8.0;
    weights[k + 3] += 3.0 / 8.0;
  }

  const std::complex<double> z(omega, gamma);
  Eigen::Matrix3cd alpha = Eigen::Matrix3cd::Zero();
  for (std::size_t n = 0; n < n_points; ++n) {
    const double t = response.time_step * static_cast<double>(n);
    alpha += (weights[n] * response.time_step * std::exp(std::complex<double>(0.0, 1.0) * z * t)) *
             response.response[n].cast<std::complex<double>>();
  }
  return alpha;
}

DipoleOrders dipole_orders(const RealTimePropagator& propagator, const std::function<Vec3(double)>& shape,
                           double amplitude, const RealTimeOptions& options) {
  if (!shape) throw std::invalid_argument("dipole_orders: no field profile");
  if (!(amplitude > 0.0)) throw std::invalid_argument("dipole_orders: amplitude must be positive");
  DipoleOrders result;
  result.time_step = options.time_step;

  auto run = [&](double lambda) {
    RealTimeOptions scaled = options;
    scaled.field = [&shape, lambda](double time) {
      const Vec3 f = shape(time);
      return Vec3{lambda * f[0], lambda * f[1], lambda * f[2]};
    };
    const RealTimeTrajectory t = propagator.propagate({0.0, 0.0, 0.0}, scaled);
    result.fock_builds += t.fock_builds;
    result.steps_converged = result.steps_converged && t.steps_converged;
    result.orthonormality_error = std::max(result.orthonormality_error, t.orthonormality_error);
    return t.dipole;
  };
  const std::vector<Vec3> p1 = run(amplitude), m1 = run(-amplitude), p2 = run(2.0 * amplitude), m2 = run(-2.0 * amplitude);
  const Vec3 mu0 = propagator.ground_state_dipole();
  const double a = amplitude;

  const std::size_t n_points = p1.size();
  result.first.resize(n_points);
  result.second.resize(n_points);
  result.third.resize(n_points);
  for (std::size_t n = 0; n < n_points; ++n)
    for (std::size_t k = 0; k < 3; ++k) {
      const double o1 = 0.5 * (p1[n][k] - m1[n][k]), o2 = 0.5 * (p2[n][k] - m2[n][k]);
      const double e1 = 0.5 * (p1[n][k] + m1[n][k]) - mu0[k], e2 = 0.5 * (p2[n][k] + m2[n][k]) - mu0[k];
      result.first[n][k] = (8.0 * o1 - o2) / (6.0 * a);
      result.third[n][k] = (o2 - 2.0 * o1) / (6.0 * a * a * a);
      result.second[n][k] = (16.0 * e1 - e2) / (12.0 * a * a);
    }
  return result;
}

}  // namespace ellidft
