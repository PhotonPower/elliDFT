// Levenberg-Marquardt least squares: known solutions, bounds, statistics.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

#include "ellipsometry/least_squares.hpp"

using namespace ellidft::ellipsometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const double inf = std::numeric_limits<double>::infinity();

Eigen::VectorXd unbounded(Eigen::Index n, double value) { return Eigen::VectorXd::Constant(n, value); }

}  // namespace

TEST_CASE("least squares: straight line with the closed solution and covariance", "[least_squares]") {
  // y = a + b x; solution of the normal equations and covariance
  // s^2 (X^T X)^-1 with s^2 = S / (n - 2).
  const Eigen::Index n = 7;
  Eigen::VectorXd x(n), y(n);
  x << 0.0, 1.0, 2.0, 3.5, 4.0, 5.5, 7.0;
  y << 1.1, 2.9, 5.2, 7.9, 9.1, 12.2, 14.8;
  Eigen::MatrixXd X(n, 2);
  X.col(0).setOnes();
  X.col(1) = x;
  const Eigen::Vector2d exact = (X.transpose() * X).ldlt().solve(X.transpose() * y);
  const double S = (y - X * exact).squaredNorm();
  const Eigen::Matrix2d covariance = S / double(n - 2) * (X.transpose() * X).inverse();

  const auto residuals = [&](const Eigen::VectorXd& p) -> Eigen::VectorXd { return y - (p[0] * Eigen::VectorXd::Ones(n) + p[1] * x); };
  Eigen::Vector2d start(-3.0, 10.0);
  const LeastSquaresResult fit = levenberg_marquardt(residuals, start, unbounded(2, -inf), unbounded(2, inf));
  CHECK(fit.converged);
  CHECK_THAT(fit.parameters[0], WithinAbs(exact[0], 1e-9));
  CHECK_THAT(fit.parameters[1], WithinAbs(exact[1], 1e-9));
  CHECK_THAT(fit.sum_of_squares, WithinRel(S, 1e-9));
  REQUIRE(fit.covariance_valid);
  CHECK((fit.covariance - covariance).cwiseAbs().maxCoeff() < 1e-8 * covariance.cwiseAbs().maxCoeff());
  CHECK_THAT(fit.standard_errors[1], WithinRel(std::sqrt(covariance(1, 1)), 1e-7));
  CHECK_THAT(fit.correlation(0, 1), WithinRel(covariance(0, 1) / std::sqrt(covariance(0, 0) * covariance(1, 1)), 1e-7));
  CHECK_THAT(fit.correlation(0, 0), WithinRel(1.0, 1e-12));
}

TEST_CASE("least squares: Rosenbrock valley", "[least_squares]") {
  // S = 100 (y - x^2)^2 + (1 - x)^2 with the minimum S = 0 at (1, 1)
  const auto residuals = [](const Eigen::VectorXd& p) -> Eigen::VectorXd {
    return Eigen::Vector2d(10.0 * (p[1] - p[0] * p[0]), 1.0 - p[0]);
  };
  const LeastSquaresResult fit = levenberg_marquardt(residuals, Eigen::Vector2d(-1.2, 1.0), unbounded(2, -inf), unbounded(2, inf));
  CHECK(fit.converged);
  CHECK_THAT(fit.parameters[0], WithinAbs(1.0, 1e-8));
  CHECK_THAT(fit.parameters[1], WithinAbs(1.0, 1e-8));
  CHECK(fit.sum_of_squares < 1e-16);
  CHECK_FALSE(fit.covariance_valid);   // as many residuals as parameters: no estimate of the scatter
}

TEST_CASE("least squares: exponential decay with parameters of different size", "[least_squares]") {
  // y = A exp(-t / tau) + c, exact data
  const double A = 250.0, tau = 0.004, c = 3.0;
  const Eigen::Index n = 40;
  Eigen::VectorXd t(n), y(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    t[i] = 0.0005 * double(i);
    y[i] = A * std::exp(-t[i] / tau) + c;
  }
  const auto residuals = [&](const Eigen::VectorXd& p) -> Eigen::VectorXd {
    Eigen::VectorXd r(n);
    for (Eigen::Index i = 0; i < n; ++i) r[i] = y[i] - (p[0] * std::exp(-t[i] / p[1]) + p[2]);
    return r;
  };
  Eigen::Vector3d lower(0.0, 1e-5, -inf), upper(inf, 1.0, inf);
  const LeastSquaresResult fit = levenberg_marquardt(residuals, Eigen::Vector3d(100.0, 0.01, 0.0), lower, upper);
  CHECK(fit.converged);
  CHECK_THAT(fit.parameters[0], WithinRel(A, 1e-8));
  CHECK_THAT(fit.parameters[1], WithinRel(tau, 1e-8));
  CHECK_THAT(fit.parameters[2], WithinRel(c, 1e-7));
}

TEST_CASE("least squares: bounds", "[least_squares]") {
  // S = (p0 - 3)^2 + (p1 + 1)^2 + ...; with p0 <= 2 and p1 >= 0 the minimum
  // of the box is at (2, 0).
  const auto residuals = [](const Eigen::VectorXd& p) -> Eigen::VectorXd {
    return Eigen::Vector3d(p[0] - 3.0, p[1] + 1.0, 0.5 * (p[0] - 3.0));
  };
  const LeastSquaresResult fit = levenberg_marquardt(residuals, Eigen::Vector2d(0.5, 4.0), Eigen::Vector2d(-inf, 0.0), Eigen::Vector2d(2.0, inf));
  CHECK(fit.converged);
  CHECK_THAT(fit.parameters[0], WithinAbs(2.0, 1e-12));
  CHECK_THAT(fit.parameters[1], WithinAbs(0.0, 1e-12));
  CHECK_THAT(fit.sum_of_squares, WithinRel(1.0 + 1.0 + 0.25, 1e-12));

  // start values outside the box are moved into it
  const LeastSquaresResult moved = levenberg_marquardt(residuals, Eigen::Vector2d(50.0, -50.0), Eigen::Vector2d(-inf, 0.0), Eigen::Vector2d(2.0, inf));
  CHECK_THAT(moved.parameters[0], WithinAbs(2.0, 1e-12));
  CHECK_THAT(moved.parameters[1], WithinAbs(0.0, 1e-12));
}

TEST_CASE("least squares: dependent parameters are reported", "[least_squares]") {
  // only the sum p0 + p1 enters
  const auto residuals = [](const Eigen::VectorXd& p) -> Eigen::VectorXd {
    return Eigen::Vector3d(p[0] + p[1] - 1.0, 2.0 * (p[0] + p[1]) - 2.1, p[0] + p[1] - 0.9);
  };
  const LeastSquaresResult fit = levenberg_marquardt(residuals, Eigen::Vector2d(0.0, 0.0), unbounded(2, -inf), unbounded(2, inf));
  CHECK_THAT(fit.parameters[0] + fit.parameters[1], WithinAbs(6.1 / 6.0, 1e-9));   // (1 + 4.2 + 0.9) / 6
  CHECK_FALSE(fit.covariance_valid);
  CHECK(std::isnan(fit.standard_errors[0]));
}

TEST_CASE("least squares: invalid input", "[least_squares]") {
  const auto one = [](const Eigen::VectorXd& p) -> Eigen::VectorXd { return Eigen::VectorXd::Constant(1, p[0]); };
  CHECK_THROWS_AS(levenberg_marquardt(one, Eigen::Vector2d(1.0, 1.0), unbounded(2, -inf), unbounded(2, inf)), std::invalid_argument);
  CHECK_THROWS_AS(levenberg_marquardt(one, unbounded(1, 1.0), unbounded(1, 2.0), unbounded(1, 1.0)), std::invalid_argument);
  CHECK_THROWS_AS(levenberg_marquardt(one, unbounded(1, 1.0), unbounded(2, 0.0), unbounded(1, 1.0)), std::invalid_argument);
  const auto bad = [](const Eigen::VectorXd&) -> Eigen::VectorXd { return Eigen::VectorXd::Constant(2, std::nan("")); };
  CHECK_THROWS_AS(levenberg_marquardt(bad, unbounded(1, 1.0), unbounded(1, -inf), unbounded(1, inf)), std::runtime_error);
}
