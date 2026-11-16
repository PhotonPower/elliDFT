#include "ellipsometry/dispersion_fit.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#include <Eigen/Dense>

#include "ellipsometry/least_squares.hpp"

namespace ellidft::ellipsometry {

namespace {

void check(const std::vector<double>& wavelengths_nm, const std::vector<double>& n, std::size_t minimum) {
  if (wavelengths_nm.size() != n.size()) throw std::invalid_argument("dispersion fit: sizes differ");
  if (n.size() < minimum) throw std::invalid_argument("dispersion fit: too few points");
  for (double w : wavelengths_nm)
    if (!(w > 0.0)) throw std::invalid_argument("dispersion fit: wavelengths must be positive");
}

}  // namespace

CauchyFit fit_cauchy(const std::vector<double>& wavelengths_nm, const std::vector<double>& n) {
  check(wavelengths_nm, n, 2);
  const auto m = static_cast<Eigen::Index>(n.size());
  const Eigen::Index terms = m >= 4 ? 3 : 2;
  Eigen::MatrixXd X(m, terms);
  Eigen::VectorXd y(m);
  for (Eigen::Index i = 0; i < m; ++i) {
    const double l2 = std::pow(wavelengths_nm[static_cast<std::size_t>(i)] / 1000.0, 2);
    X(i, 0) = 1.0;
    X(i, 1) = 1.0 / l2;
    if (terms == 3) X(i, 2) = 1.0 / (l2 * l2);
    y[i] = n[static_cast<std::size_t>(i)];
  }
  const Eigen::VectorXd c = X.colPivHouseholderQr().solve(y);
  CauchyFit fit;
  fit.A = c[0];
  fit.B = c[1];
  fit.C = terms == 3 ? c[2] : 0.0;
  const Eigen::VectorXd r = y - X * c;
  fit.rms = std::sqrt(r.squaredNorm() / static_cast<double>(m));
  fit.largest = r.cwiseAbs().maxCoeff();
  return fit;
}

SellmeierFit fit_sellmeier(const std::vector<double>& wavelengths_nm, const std::vector<double>& n) {
  check(wavelengths_nm, n, 3);
  const auto m = static_cast<Eigen::Index>(n.size());
  Eigen::MatrixXd X(m, 2);
  Eigen::VectorXd y(m), l2(m), data(m);
  for (Eigen::Index i = 0; i < m; ++i) {
    const double value = n[static_cast<std::size_t>(i)];
    if (!(value > 1.0)) throw std::invalid_argument("dispersion fit: the Sellmeier model needs n > 1");
    l2[i] = std::pow(wavelengths_nm[static_cast<std::size_t>(i)] / 1000.0, 2);
    data[i] = value;
    X(i, 0) = 1.0;
    X(i, 1) = -1.0 / l2[i];
    y[i] = 1.0 / (value * value - 1.0);
  }
  const Eigen::Vector2d line = X.colPivHouseholderQr().solve(y);   // 1/B and C/B
  Eigen::Vector2d start(1.0 / line[0], line[1] / line[0]);
  const double smallest_l2 = l2.minCoeff();
  if (!(start[0] > 0.0)) start[0] = 1.0;
  if (!(start[1] >= 0.0) || start[1] >= smallest_l2) start[1] = 0.0;

  const ResidualFunction residuals = [&](const Eigen::VectorXd& p) -> Eigen::VectorXd {
    Eigen::VectorXd r(m);
    for (Eigen::Index i = 0; i < m; ++i) r[i] = data[i] - std::sqrt(1.0 + p[0] * l2[i] / (l2[i] - p[1]));
    return r;
  };
  // the pole must stay below the shortest wavelength of the data
  const Eigen::Vector2d lower(0.0, 0.0), upper(std::numeric_limits<double>::infinity(), 0.999 * smallest_l2);
  const LeastSquaresResult result = levenberg_marquardt(residuals, start, lower, upper);
  SellmeierFit fit;
  fit.B = result.parameters[0];
  fit.C_um2 = result.parameters[1];
  fit.converged = result.converged;
  const Eigen::VectorXd r = residuals(result.parameters);
  fit.rms = std::sqrt(r.squaredNorm() / static_cast<double>(m));
  fit.largest = r.cwiseAbs().maxCoeff();
  return fit;
}

}  // namespace ellidft::ellipsometry
