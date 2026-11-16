#pragma once

// Non-linear least squares with the Levenberg-Marquardt method, with simple
// bounds on the parameters.

#include <functional>
#include <string>

#include <Eigen/Dense>

namespace ellidft::ellipsometry {

// Residuals r_i(p); minimised is S(p) = sum_i r_i(p)^2.
using ResidualFunction = std::function<Eigen::VectorXd(const Eigen::VectorXd&)>;

struct LeastSquaresOptions {
  int max_iterations = 200;
  double reduction_tolerance = 1e-12;  // stop when S decreases by less than this, relative
  double step_tolerance = 1e-10;       // stop when the step is smaller than this, relative
  double gradient_tolerance = 1e-12;   // stop when the scaled gradient is smaller than this
  double initial_damping = 1e-3;
  double derivative_step = 1e-6;       // relative step of the central differences
};

struct LeastSquaresResult {
  Eigen::VectorXd parameters;
  double sum_of_squares = 0.0;
  int iterations = 0;
  int evaluations = 0;
  bool converged = false;
  std::string message;

  // Statistics at the solution. With n residuals and m parameters:
  //   covariance = S / (n - m) * (J^T J)^-1,
  // i.e. the residuals are taken to have a common, unknown standard
  // deviation that is estimated from the fit itself. standard_errors are the
  // square roots of its diagonal, correlation its normalised form. If J^T J is
  // singular (parameters not independent, or a parameter without influence)
  // `covariance_valid` is false and the three fields hold NaN.
  bool covariance_valid = false;
  Eigen::MatrixXd covariance, correlation;
  Eigen::VectorXd standard_errors;
};

// `lower` and `upper` may hold -inf / +inf. The start values are moved into
// the bounds if necessary.
LeastSquaresResult levenberg_marquardt(const ResidualFunction& residuals, Eigen::VectorXd start,
                                       const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
                                       const LeastSquaresOptions& options = {});

}  // namespace ellidft::ellipsometry
