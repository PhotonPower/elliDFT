#include "scf/diis.hpp"

#include <cmath>
#include <stdexcept>

namespace ellidft {

Diis::Diis(std::size_t max_vectors) : max_vectors_(max_vectors) {
  if (max_vectors == 0) throw std::invalid_argument("Diis: at least one vector is needed");
}

Eigen::MatrixXd Diis::extrapolate(const Eigen::MatrixXd& fock, const Eigen::MatrixXd& error) {
  focks_.push_back(fock);
  errors_.push_back(error);
  if (focks_.size() > max_vectors_) {
    focks_.pop_front();
    errors_.pop_front();
  }

  // If the linear system becomes numerically singular (nearly parallel error
  // matrices), the oldest pairs are discarded and the system is solved again.
  while (true) {
    const auto m = static_cast<Eigen::Index>(focks_.size());
    if (m == 1) return focks_.back();

    // Hernández Vera 2025, Eq. 69: B_ij = <e_i, e_j> (Frobenius product),
    //   (B 1; 1^T 0) (c; lambda) = (0; 1).
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(m + 1, m + 1);
    for (Eigen::Index i = 0; i < m; ++i)
      for (Eigen::Index j = 0; j <= i; ++j)
        A(i, j) = A(j, i) = errors_[static_cast<std::size_t>(i)].cwiseProduct(errors_[static_cast<std::size_t>(j)]).sum();
    // Scaling B changes lambda only; it keeps the system well scaled when
    // the errors become tiny near convergence.
    const double scale = A.topLeftCorner(m, m).diagonal().maxCoeff();
    if (scale > 0.0) A.topLeftCorner(m, m) /= scale;
    A.block(0, m, m, 1).setOnes();
    A.block(m, 0, 1, m).setOnes();
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m + 1);
    rhs(m) = 1.0;

    const Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    const Eigen::VectorXd solution = lu.solve(rhs);
    const bool usable = lu.isInvertible() && solution.allFinite() &&
                        (A * solution - rhs).cwiseAbs().maxCoeff() < 1e-8 &&
                        solution.head(m).cwiseAbs().maxCoeff() < 1e6;
    if (!usable) {
      focks_.pop_front();
      errors_.pop_front();
      continue;
    }

    // Eq. 66: F = sum_i c_i F_i
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(fock.rows(), fock.cols());
    for (Eigen::Index i = 0; i < m; ++i) result += solution(i) * focks_[static_cast<std::size_t>(i)];
    return result;
  }
}

}  // namespace ellidft
