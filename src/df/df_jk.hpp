#pragma once

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "scf/jk.hpp"

namespace ellidft {

// Coulomb and exchange matrices with density fitting in the Coulomb metric
// (RI-J and RI-K):
//
//   (mu nu|lambda sigma) ~ sum_PQ (mu nu|P) [V^-1]_PQ (Q|lambda sigma),  V_PQ = (P|Q).
//
// With the Cholesky factor V = L L^T and B^Q_mu,nu = sum_P [L^-1]_QP (P|mu nu):
//
//   J_mu,nu = sum_Q B^Q_mu,nu  sum_{lambda sigma} B^Q_lambda,sigma D_lambda,sigma
//   K_mu,nu = sum_Q sum_{lambda sigma} B^Q_mu,lambda D_lambda,sigma B^Q_nu,sigma
//
// The fitted three-index quantity B is held in memory, packed over
// mu >= nu: n(n+1)/2 * n_aux numbers.
class DensityFittedJK final : public JKBuilder {
 public:
  // Throws std::runtime_error if the metric is not positive definite
  // (linearly dependent auxiliary basis).
  DensityFittedJK(const BasisSet& basis, const BasisSet& auxiliary);

  void build(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
             std::vector<Eigen::MatrixXd>* K) override;

  // Long-range exchange fitted in the metric of the attenuated operator:
  // both the three-centre integrals and the metric use erf(omega r)/r, as in
  // the reference code. The fitted quantity is built on the first call for a
  // given omega and kept.
  void build_long_range_exchange(const std::vector<Eigen::MatrixXd>& densities, double omega,
                                 std::vector<Eigen::MatrixXd>& K_long_range) override;

  Eigen::Index n_auxiliary() const { return B_.cols(); }

  // For protocols and for the decision whether the three-index quantity
  // still fits into memory: bytes held (both operators, once the long-range
  // one has been built) and wall time spent building them.
  double memory_bytes() const { return 8.0 * static_cast<double>(B_.size() + B_long_range_.size()); }
  double build_seconds() const { return build_seconds_; }

 private:
  static Eigen::MatrixXd fitted(const BasisSet& basis, const BasisSet& auxiliary, double omega);
  void exchange(const Eigen::MatrixXd& B, const std::vector<Eigen::MatrixXd>& densities,
                std::vector<Eigen::MatrixXd>& K) const;

  Eigen::Index n_;
  BasisSet basis_, auxiliary_;
  Eigen::MatrixXd B_;  // B_(mu (mu+1)/2 + nu, Q)
  double long_range_omega_ = 0.0;
  Eigen::MatrixXd B_long_range_;
  double build_seconds_ = 0.0;
};

}  // namespace ellidft
