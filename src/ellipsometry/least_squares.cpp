#include "ellipsometry/least_squares.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ellidft::ellipsometry {

namespace {

// Jacobian by central differences; one-sided where a bound is in the way.
Eigen::MatrixXd jacobian(const ResidualFunction& residuals, const Eigen::VectorXd& p, const Eigen::VectorXd& r0,
                         const Eigen::VectorXd& lower, const Eigen::VectorXd& upper, double relative_step,
                         int& evaluations) {
  Eigen::MatrixXd J(r0.size(), p.size());
  for (Eigen::Index j = 0; j < p.size(); ++j) {
    const double h = relative_step * std::max(std::abs(p[j]), 1e-3);
    Eigen::VectorXd plus = p, minus = p;
    plus[j] = std::min(p[j] + h, upper[j]);
    minus[j] = std::max(p[j] - h, lower[j]);
    const double width = plus[j] - minus[j];
    if (!(width > 0.0)) {  // parameter pinned by coinciding bounds
      J.col(j).setZero();
      continue;
    }
    const Eigen::VectorXd r_plus = plus[j] == p[j] ? r0 : residuals(plus);
    const Eigen::VectorXd r_minus = minus[j] == p[j] ? r0 : residuals(minus);
    evaluations += (plus[j] != p[j]) + (minus[j] != p[j]);
    J.col(j) = (r_plus - r_minus) / width;
  }
  return J;
}

}  // namespace

LeastSquaresResult levenberg_marquardt(const ResidualFunction& residuals, Eigen::VectorXd start,
                                       const Eigen::VectorXd& lower, const Eigen::VectorXd& upper,
                                       const LeastSquaresOptions& options) {
  const Eigen::Index m = start.size();
  if (lower.size() != m || upper.size() != m) throw std::invalid_argument("levenberg_marquardt: sizes of the bounds");
  for (Eigen::Index j = 0; j < m; ++j) {
    if (lower[j] > upper[j]) throw std::invalid_argument("levenberg_marquardt: lower bound above upper bound");
    start[j] = std::clamp(start[j], lower[j], upper[j]);
  }

  LeastSquaresResult result;
  Eigen::VectorXd p = start;
  Eigen::VectorXd r = residuals(p);
  result.evaluations = 1;
  const Eigen::Index n = r.size();
  if (n < m) throw std::invalid_argument("levenberg_marquardt: fewer residuals than parameters");
  if (!r.allFinite()) throw std::runtime_error("levenberg_marquardt: residuals are not finite at the start values");
  double S = r.squaredNorm();
  double damping = options.initial_damping;
  Eigen::MatrixXd J;

  // Levenberg-Marquardt in the form with the diagonal of J^T J (Wikipedia,
  // "Levenberg-Marquardt algorithm", section "The solution"):
  //   [J^T J + lambda diag(J^T J)] delta = J^T [y - f(p)] = -J^T r.
  // Damping: decrease by a factor 3 after a successful step, increase by a
  // factor 2 after a rejected one (same page, "Choice of damping parameter").
  // Bounds: the trial point is projected onto the box; components that sit
  // on a bound and would leave it are taken out of the linear system.
  for (result.iterations = 0; result.iterations < options.max_iterations;) {
    J = jacobian(residuals, p, r, lower, upper, options.derivative_step, result.evaluations);
    const Eigen::VectorXd g = J.transpose() * r;
    Eigen::MatrixXd A = J.transpose() * J;

    std::vector<Eigen::Index> free;
    for (Eigen::Index j = 0; j < m; ++j) {
      const bool at_lower = p[j] <= lower[j] && g[j] > 0.0;   // step -g would leave the box
      const bool at_upper = p[j] >= upper[j] && g[j] < 0.0;
      if (!at_lower && !at_upper && A(j, j) > 0.0) free.push_back(j);
    }
    // Scaled gradient: cosine of the angle between the residual vector and
    // the columns of J; zero at a minimum.
    double scaled_gradient = 0.0;
    for (Eigen::Index j : free)
      scaled_gradient = std::max(scaled_gradient, std::abs(g[j]) / (std::sqrt(A(j, j)) * std::max(std::sqrt(S), 1e-300)));
    if (free.empty() || scaled_gradient <= options.gradient_tolerance || S == 0.0) {
      result.converged = true;
      result.message = "gradient vanishes";
      break;
    }

    const Eigen::Index nf = static_cast<Eigen::Index>(free.size());
    Eigen::MatrixXd Af(nf, nf);
    Eigen::VectorXd gf(nf);
    for (Eigen::Index a = 0; a < nf; ++a) {
      gf[a] = g[free[static_cast<std::size_t>(a)]];
      for (Eigen::Index b = 0; b < nf; ++b)
        Af(a, b) = A(free[static_cast<std::size_t>(a)], free[static_cast<std::size_t>(b)]);
    }

    bool accepted = false;
    while (!accepted) {
      Eigen::MatrixXd damped = Af;
      damped.diagonal() += damping * Af.diagonal();
      const Eigen::VectorXd delta = damped.ldlt().solve(-gf);
      Eigen::VectorXd trial = p;
      for (Eigen::Index a = 0; a < nf; ++a) {
        const Eigen::Index j = free[static_cast<std::size_t>(a)];
        trial[j] = std::clamp(p[j] + delta[a], lower[j], upper[j]);
      }
      const double step = (trial - p).norm();
      if (!delta.allFinite() || step <= options.step_tolerance * (p.norm() + options.step_tolerance)) {
        result.converged = true;
        result.message = "step vanishes";
        break;
      }
      const Eigen::VectorXd r_trial = residuals(trial);
      ++result.evaluations;
      const double S_trial = r_trial.allFinite() ? r_trial.squaredNorm() : std::numeric_limits<double>::infinity();
      if (S_trial < S) {
        const double reduction = (S - S_trial) / S;
        p = trial;
        r = r_trial;
        S = S_trial;
        damping = std::max(damping / 3.0, 1e-15);
        accepted = true;
        ++result.iterations;
        if (reduction <= options.reduction_tolerance) {
          result.converged = true;
          result.message = "sum of squares no longer decreases";
        }
      } else {
        damping *= 2.0;
        if (damping > 1e15) {
          result.converged = true;
          result.message = "no further reduction possible";
          break;
        }
      }
    }
    if (result.converged) break;
  }
  if (!result.converged) result.message = "maximum number of iterations reached";

  result.parameters = p;
  result.sum_of_squares = S;

  // Covariance from the Jacobian at the solution. The matrix is brought to
  // unit diagonal first so that the test for singularity does not depend on
  // the units of the parameters.
  const double nan = std::numeric_limits<double>::quiet_NaN();
  result.covariance = Eigen::MatrixXd::Constant(m, m, nan);
  result.correlation = Eigen::MatrixXd::Constant(m, m, nan);
  result.standard_errors = Eigen::VectorXd::Constant(m, nan);
  if (n > m) {
    J = jacobian(residuals, p, r, lower, upper, options.derivative_step, result.evaluations);
    const Eigen::MatrixXd A = J.transpose() * J;
    const Eigen::VectorXd d = A.diagonal().cwiseSqrt();
    if ((d.array() > 0.0).all()) {
      const Eigen::MatrixXd scaled = d.cwiseInverse().asDiagonal() * A * d.cwiseInverse().asDiagonal();
      const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(scaled);
      if (eig.eigenvalues().minCoeff() > 1e-12) {
        const Eigen::MatrixXd inverse = eig.eigenvectors() * eig.eigenvalues().cwiseInverse().asDiagonal() *
                                        eig.eigenvectors().transpose();
        const double variance = S / static_cast<double>(n - m);
        result.covariance = variance * (d.cwiseInverse().asDiagonal() * inverse * d.cwiseInverse().asDiagonal());
        result.standard_errors = result.covariance.diagonal().cwiseSqrt();
        const Eigen::VectorXd s = inverse.diagonal().cwiseSqrt();
        result.correlation = s.cwiseInverse().asDiagonal() * inverse * s.cwiseInverse().asDiagonal();
        result.covariance_valid = true;
      }
    }
  }
  return result;
}

}  // namespace ellidft::ellipsometry
