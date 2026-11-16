#include "df/df_jk.hpp"

#include <chrono>

#include <stdexcept>

#include "core/parallel.hpp"
#include "df/df_integrals.hpp"

namespace ellidft {

namespace {
double wall_seconds() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
}  // namespace

DensityFittedJK::DensityFittedJK(const BasisSet& basis, const BasisSet& auxiliary)
    : n_(basis.n_functions()), basis_(basis), auxiliary_(auxiliary) {
  const double start = wall_seconds();
  B_ = fitted(basis, auxiliary, 0.0);
  build_seconds_ = wall_seconds() - start;
}

Eigen::MatrixXd DensityFittedJK::fitted(const BasisSet& basis, const BasisSet& auxiliary, double omega) {
  const Eigen::MatrixXd V = two_center_coulomb(auxiliary, omega);
  const Eigen::LLT<Eigen::MatrixXd> cholesky(V);
  if (cholesky.info() != Eigen::Success)
    throw std::runtime_error("DensityFittedJK: the Coulomb metric of the auxiliary basis is not positive definite");
  // B = A L^-T with A_(mu nu),P = (mu nu|P): B L^T = A, solved in place from
  // the right with the upper triangle L^T. No second copy of the three-index
  // quantity is made (the earlier form L B^T = A^T held up to three).
  Eigen::MatrixXd A = three_center_coulomb(basis, auxiliary, omega);
  cholesky.matrixU().solveInPlace<Eigen::OnTheRight>(A);
  return A;
}

void DensityFittedJK::build_long_range_exchange(const std::vector<Eigen::MatrixXd>& densities, double omega,
                                                std::vector<Eigen::MatrixXd>& K_long_range) {
  if (!(omega > 0.0)) throw std::invalid_argument("DensityFittedJK: long-range exchange needs omega > 0");
  for (const Eigen::MatrixXd& D : densities)
    if (D.rows() != n_ || D.cols() != n_)
      throw std::invalid_argument("DensityFittedJK: density matrix has the wrong dimension");
  if (B_long_range_.size() == 0 || omega != long_range_omega_) {
    const double start = wall_seconds();
    B_long_range_ = fitted(basis_, auxiliary_, omega);
    long_range_omega_ = omega;
    build_seconds_ += wall_seconds() - start;
  }
  exchange(B_long_range_, densities, K_long_range);
}

void DensityFittedJK::exchange(const Eigen::MatrixXd& B, const std::vector<Eigen::MatrixXd>& densities,
                               std::vector<Eigen::MatrixXd>& K) const {
  const Eigen::Index n_aux = B.cols();
  const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(n_, n_);
  const int n_threads = max_threads();
  std::vector<std::vector<Eigen::MatrixXd>> Kt(static_cast<std::size_t>(n_threads));
#pragma omp parallel num_threads(n_threads)
  {
    std::vector<Eigen::MatrixXd>& Kl = Kt[static_cast<std::size_t>(thread_index())];
    Kl.assign(densities.size(), zero);
    Eigen::MatrixXd BQ(n_, n_), half(n_, n_);
#pragma omp for schedule(static)
    for (Eigen::Index Q = 0; Q < n_aux; ++Q) {
      for (Eigen::Index i = 0; i < n_; ++i)
        for (Eigen::Index j = 0; j <= i; ++j) BQ(i, j) = BQ(j, i) = B(i * (i + 1) / 2 + j, Q);
      for (std::size_t d = 0; d < densities.size(); ++d) {
        half.noalias() = BQ * densities[d];
        Kl[d].noalias() += half * BQ;  // B^Q D B^Q, B^Q symmetric
      }
    }
  }
  K.assign(densities.size(), zero);
  for (std::size_t d = 0; d < densities.size(); ++d)
    for (const auto& thread_K : Kt)
      if (!thread_K.empty()) K[d] += thread_K[d];
}

void DensityFittedJK::build(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
                            std::vector<Eigen::MatrixXd>* K) {
  for (const Eigen::MatrixXd& D : densities)
    if (D.rows() != n_ || D.cols() != n_)
      throw std::invalid_argument("DensityFittedJK: density matrix has the wrong dimension");
  const Eigen::Index n_pairs = B_.rows();

  if (J != nullptr) {
    J->clear();
    for (const Eigen::MatrixXd& D : densities) {
      // gamma_Q = sum_{mu nu} B^Q_mu,nu D_mu,nu; the packed sum runs over
      // mu >= nu with the two off-diagonal elements of D added up.
      Eigen::VectorXd packed(n_pairs);
      for (Eigen::Index i = 0; i < n_; ++i)
        for (Eigen::Index j = 0; j <= i; ++j) packed(i * (i + 1) / 2 + j) = (i == j) ? D(i, i) : D(i, j) + D(j, i);
      const Eigen::VectorXd gamma = B_.transpose() * packed;
      const Eigen::VectorXd j_packed = B_ * gamma;
      Eigen::MatrixXd Jm(n_, n_);
      for (Eigen::Index i = 0; i < n_; ++i)
        for (Eigen::Index j = 0; j <= i; ++j) Jm(i, j) = Jm(j, i) = j_packed(i * (i + 1) / 2 + j);
      J->push_back(std::move(Jm));
    }
  }

  if (K != nullptr) exchange(B_, densities, *K);
}

}  // namespace ellidft
