#include "response/polarizability.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace ellidft {

namespace {

using Complex = std::complex<double>;

// Appends the columns of `candidates` to the orthonormal basis `basis`,
// after projecting out the basis (twice, for numerical orthogonality) and
// dropping what is left with a norm below `threshold`.
Eigen::Index extend_basis(Eigen::MatrixXd& basis, const Eigen::MatrixXd& candidates, double threshold = 1e-8) {
  Eigen::Index added = 0;
  for (Eigen::Index c = 0; c < candidates.cols(); ++c) {
    Eigen::VectorXd v = candidates.col(c);
    const double initial = v.norm();
    if (!(initial > 0.0)) continue;
    v /= initial;
    for (int pass = 0; pass < 2; ++pass)
      if (basis.cols() > 0) v -= basis * (basis.transpose() * v);
    const double norm = v.norm();
    if (norm < threshold) continue;
    basis.conservativeResize(candidates.rows(), basis.cols() + 1);
    basis.col(basis.cols() - 1) = v / norm;
    ++added;
  }
  return added;
}

}  // namespace

Polarizability polarizability(const ResponseOperator& op, double omega, double gamma,
                              const LinearResponseOptions& options) {
  Polarizability result;
  result.omega = omega;
  result.gamma = gamma;

  const Eigen::Index N = op.size();
  const Eigen::MatrixXd mu = op.dipole_vectors();          // N x 3
  const Eigen::VectorXd& d = op.orbital_energy_differences();
  const Complex z(omega, gamma);
  const bool is_static = omega == 0.0 && gamma == 0.0;
  const bool is_complex = gamma != 0.0;

  // Preconditioner: the equations with A+B and A-B replaced by the orbital
  // energy differences decouple into 2x2 systems per orbital pair,
  //   (d  -z; -z  d) (t1; t2) = (r1; r2)  =>  t1 = (d r1 + z r2) / (d^2 - z^2), ...
  auto precondition = [&](const Eigen::VectorXcd& r1, const Eigen::VectorXcd& r2, Eigen::VectorXcd& t1,
                          Eigen::VectorXcd& t2) {
    t1.resize(N);
    t2.resize(N);
    for (Eigen::Index k = 0; k < N; ++k) {
      Complex denominator = d(k) * d(k) - z * z;
      if (std::abs(denominator) < 1e-8) denominator = Complex(1e-8, 0.0);  // near a pole of the model
      t1(k) = (d(k) * r1(k) + z * r2(k)) / denominator;
      t2(k) = (d(k) * r2(k) + z * r1(k)) / denominator;
    }
  };
  auto real_and_imaginary_parts = [&](const std::vector<Eigen::VectorXcd>& vectors) {
    Eigen::MatrixXd out(N, static_cast<Eigen::Index>(vectors.size()) * (is_complex ? 2 : 1));
    Eigen::Index c = 0;
    for (const Eigen::VectorXcd& v : vectors) {
      out.col(c++) = v.real();
      if (is_complex) out.col(c++) = v.imag();
    }
    return out;
  };

  // Start vectors: the preconditioned right-hand sides (2 mu, 0).
  Eigen::MatrixXd basis(N, 0);
  {
    std::vector<Eigen::VectorXcd> start;
    for (Eigen::Index l = 0; l < 3; ++l) {
      Eigen::VectorXcd t1, t2;
      precondition((2.0 * mu.col(l)).cast<Complex>(), Eigen::VectorXcd::Zero(N), t1, t2);
      start.push_back(t1);
      if (!is_static) start.push_back(t2);
    }
    extend_basis(basis, real_and_imaginary_parts(start));
  }
  if (basis.cols() == 0) {  // no dipole-allowed response at all
    result.converged = true;
    return result;
  }

  Eigen::MatrixXd plus(N, 0), minus(N, 0);  // (A+B) b and (A-B) b for the basis vectors
  Eigen::Index processed = 0;
  Eigen::MatrixXcd Z(N, 3), W(N, 3);

  for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
    result.iterations = iteration;
    // Products for the new basis vectors.
    Eigen::MatrixXd new_plus, new_minus;
    op.apply(basis.rightCols(basis.cols() - processed), &new_plus, &new_minus);
    plus.conservativeResize(N, basis.cols());
    minus.conservativeResize(N, basis.cols());
    plus.rightCols(new_plus.cols()) = new_plus;
    minus.rightCols(new_minus.cols()) = new_minus;
    processed = basis.cols();

    // Reduced equations: (Sp  -z; -z  Sm) (zeta; w) = (2 b^T mu; 0).
    const Eigen::Index m = basis.cols();
    const Eigen::MatrixXd Sp = basis.transpose() * plus, Sm = basis.transpose() * minus;
    Eigen::MatrixXcd reduced = Eigen::MatrixXcd::Zero(2 * m, 2 * m);
    reduced.topLeftCorner(m, m) = Sp.cast<Complex>();
    reduced.bottomRightCorner(m, m) = Sm.cast<Complex>();
    reduced.topRightCorner(m, m) = -z * Eigen::MatrixXcd::Identity(m, m);
    reduced.bottomLeftCorner(m, m) = -z * Eigen::MatrixXcd::Identity(m, m);
    Eigen::MatrixXcd rhs = Eigen::MatrixXcd::Zero(2 * m, 3);
    rhs.topRows(m) = (2.0 * basis.transpose() * mu).cast<Complex>();
    const Eigen::MatrixXcd solution = reduced.fullPivLu().solve(rhs);
    const Eigen::MatrixXcd zeta = solution.topRows(m), w = solution.bottomRows(m);
    Z = basis.cast<Complex>() * zeta;
    W = basis.cast<Complex>() * w;

    // Residuals of the full equations.
    const Eigen::MatrixXcd R1 = plus.cast<Complex>() * zeta - z * W - (2.0 * mu).cast<Complex>();
    const Eigen::MatrixXcd R2 = minus.cast<Complex>() * w - z * Z;
    std::vector<Eigen::VectorXcd> corrections;
    result.residual = 0.0;
    for (Eigen::Index l = 0; l < 3; ++l) {
      const double scale = (2.0 * mu.col(l)).norm();
      if (!(scale > 0.0)) continue;
      const double residual = std::max(R1.col(l).norm(), R2.col(l).norm()) / scale;
      result.residual = std::max(result.residual, residual);
      if (residual < options.tolerance) continue;
      Eigen::VectorXcd t1, t2;
      precondition(R1.col(l), R2.col(l), t1, t2);
      corrections.push_back(t1);
      if (!is_static) corrections.push_back(t2);
    }
    if (corrections.empty()) {
      result.converged = true;
      break;
    }
    if (extend_basis(basis, real_and_imaginary_parts(corrections)) == 0) break;  // subspace exhausted
    if (basis.cols() > N) break;
  }

  // alpha_kl = 2 mu^k . Z^l
  const Eigen::Matrix3cd alpha = 2.0 * mu.transpose().cast<Complex>() * Z;
  result.real = alpha.real();
  result.imag = alpha.imag();
  return result;
}

}  // namespace ellidft
