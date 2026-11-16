#pragma once

// From real-time propagations to the polarisability.

#include <array>
#include <complex>
#include <functional>
#include <vector>

#include <Eigen/Dense>

#include "realtime/propagator.hpp"

namespace ellidft {

struct KickOptions {
  double strength = 1e-4;   // kappa, atomic units (field strength times time)
  // Two propagations per direction with +kappa and -kappa; the response is
  // their difference divided by 2 kappa. This removes all even orders in the
  // field (hyperpolarisability) and any residual motion of an imperfectly
  // converged ground state. With false: one propagation, difference to the
  // ground-state dipole divided by kappa.
  bool two_sided = true;
  std::array<bool, 3> directions{true, true, true};   // pulses along x, y, z
};

// Linear response of the dipole moment to a delta pulse:
//   response[n](i, j) = d mu_i(t_n) / d kappa_j,   t_n = n * time_step.
// Columns of directions that were not requested are zero.
struct DipoleResponse {
  double time_step = 0.0;
  std::vector<Eigen::Matrix3d> response;

  // Diagnostics over all propagations.
  int propagations = 0;
  int fock_builds = 0;
  bool steps_converged = true;
  double energy_drift = 0.0;
  double orthonormality_error = 0.0;

  double total_time() const { return time_step * static_cast<double>(response.empty() ? 0 : response.size() - 1); }
};

DipoleResponse dipole_response(const RealTimePropagator& propagator, const RealTimeOptions& options = {},
                               const KickOptions& kick = {});

// Polarisability from the response function. For a pulse F(t) = kappa
// delta(t) the Fourier transform of the field is the constant kappa, so
//
//   alpha_ij(omega + i gamma) = integral_0^T  [d mu_i(t) / d kappa_j]  e^{i (omega + i gamma) t} dt
//
// (damped Fourier, i.e. Laplace, transform: Hauge et al. 2023, Eq. 8;
// alpha = d mu(omega) / d E(omega): Herbert 2022, Eq. 82). gamma > 0 is the
// same damping as in polarizability() of the module response. The integral
// is evaluated with Simpson's rule on the time grid and ends at the last
// time T; the neglected tail is of relative size exp(-gamma T).
Eigen::Matrix3cd polarizability(const DipoleResponse& response, double omega, double gamma);

// Dipole moment order by order in the strength of a field with a given time
// profile. For the field lambda * shape(t) the dipole moment is expanded as
//
//   mu(t; lambda) = mu0 + lambda mu1(t) + lambda^2 mu2(t) + lambda^3 mu3(t) + ...
//
// mu1 is the linear response to the profile, mu2 and mu3 the responses of
// second and third order (in the static limit alpha.u, beta:uu / 2 and
// gamma:uuu / 6). They are obtained from four propagations with
// lambda = +-amplitude and +-2 amplitude by finite differences: with the odd
// and even parts o(l) = [mu(l) - mu(-l)] / 2, e(l) = [mu(l) + mu(-l)] / 2 - mu0,
//
//   mu1 = [8 o(a) - o(2a)] / (6 a),      mu3 = [o(2a) - 2 o(a)] / (6 a^3),
//   mu2 = [16 e(a) - e(2a)] / (12 a^2).
//
// The remaining errors are of the order a^4 mu5 (in mu1), a^2 mu5 (in mu3)
// and a^4 mu6 (in mu2): the amplitude must be small enough for the expansion
// and large enough for mu3 to stand out of the numerical noise.
struct DipoleOrders {
  double time_step = 0.0;
  std::vector<Vec3> first, second, third;   // mu1, mu2, mu3 at t_n = n * time_step

  int fock_builds = 0;
  bool steps_converged = true;
  double orthonormality_error = 0.0;
};

// `shape` is the time profile of the field (vector, atomic units of field
// strength for lambda = 1); options.field is ignored.
DipoleOrders dipole_orders(const RealTimePropagator& propagator, const std::function<Vec3(double)>& shape,
                           double amplitude, const RealTimeOptions& options = {});

}  // namespace ellidft
