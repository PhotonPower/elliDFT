#include "response/response_operator.hpp"

#include <stdexcept>
#include <vector>

#include "integrals/one_electron.hpp"
#include "scf/xc_matrix.hpp"

namespace ellidft {

ResponseOperator::ResponseOperator(const Molecule& molecule, const BasisSet& basis, const Method& method,
                                   const ScfResult& scf, JKBuilder& jk)
    : basis_(basis), functional_(method.functional), jk_(jk), density_(scf.density) {
  if (!scf.converged) throw std::invalid_argument("ResponseOperator: the SCF is not converged");
  n_occ_ = scf.n_occupied;
  n_virt_ = scf.coefficients.cols() - n_occ_;
  if (n_occ_ <= 0 || n_virt_ <= 0)
    throw std::invalid_argument("ResponseOperator: needs occupied and virtual orbitals");
  C_occ_ = scf.coefficients.leftCols(n_occ_);
  C_virt_ = scf.coefficients.rightCols(n_virt_);
  differences_.resize(n_occ_ * n_virt_);
  for (Eigen::Index i = 0; i < n_occ_; ++i)
    for (Eigen::Index a = 0; a < n_virt_; ++a)
      differences_(i * n_virt_ + a) = scf.orbital_energies(n_occ_ + a) - scf.orbital_energies(i);
  if (functional_.has_density_part()) {
    xc_ = make_libxc_evaluator(functional_);
    grid_ = build_molecular_grid(molecule, method.grid);
  }
}

void ResponseOperator::apply(const Eigen::MatrixXd& V, Eigen::MatrixXd* plus, Eigen::MatrixXd* minus) const {
  if (V.rows() != size()) throw std::invalid_argument("ResponseOperator::apply: vector has the wrong length");
  ++applications_;
  const Eigen::Index n_vectors = V.cols();
  const double c_x = functional_.exact_exchange, c_lr = functional_.long_range_exchange;

  // Transition densities P = C_occ V C_virt^T (not symmetric).
  std::vector<Eigen::MatrixXd> P;
  P.reserve(static_cast<std::size_t>(n_vectors));
  for (Eigen::Index k = 0; k < n_vectors; ++k) {
    const Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> Vk(
        V.col(k).data(), n_occ_, n_virt_);
    P.push_back(C_occ_ * Vk * C_virt_.transpose());
  }

  std::vector<Eigen::MatrixXd> J, K, K_lr, K_xc;
  const bool need_J = plus != nullptr;
  if (need_J || c_x != 0.0) jk_.build(P, need_J ? &J : nullptr, c_x != 0.0 ? &K : nullptr);
  if (c_lr != 0.0) jk_.build_long_range_exchange(P, functional_.omega, K_lr);
  if (plus != nullptr && xc_) K_xc = integrate_xc_kernel(basis_, grid_, *xc_, density_, P);

  auto to_vector = [&](const Eigen::MatrixXd& ao, Eigen::Ref<Eigen::VectorXd> out) {
    const Eigen::MatrixXd mo = C_occ_.transpose() * ao * C_virt_;
    for (Eigen::Index i = 0; i < n_occ_; ++i)
      for (Eigen::Index a = 0; a < n_virt_; ++a) out(i * n_virt_ + a) = mo(i, a);
  };

  if (plus != nullptr) plus->resize(size(), n_vectors);
  if (minus != nullptr) minus->resize(size(), n_vectors);
  Eigen::VectorXd contribution(size());
  for (Eigen::Index k = 0; k < n_vectors; ++k) {
    const auto i = static_cast<std::size_t>(k);
    if (plus != nullptr) {
      Eigen::MatrixXd G = 4.0 * J[i];
      if (xc_) G += 4.0 * K_xc[i];
      if (c_x != 0.0) G -= c_x * (K[i] + K[i].transpose());
      if (c_lr != 0.0) G -= c_lr * (K_lr[i] + K_lr[i].transpose());
      to_vector(G, contribution);
      plus->col(k) = differences_.cwiseProduct(V.col(k)) + contribution;
    }
    if (minus != nullptr) {
      minus->col(k) = differences_.cwiseProduct(V.col(k));
      if (c_x != 0.0 || c_lr != 0.0) {
        Eigen::MatrixXd G = Eigen::MatrixXd::Zero(density_.rows(), density_.cols());
        if (c_x != 0.0) G -= c_x * (K[i] - K[i].transpose());
        if (c_lr != 0.0) G -= c_lr * (K_lr[i] - K_lr[i].transpose());
        to_vector(G, contribution);
        minus->col(k) += contribution;
      }
    }
  }
}

Eigen::MatrixXd ResponseOperator::dipole_vectors() const {
  const auto dipole = dipole_matrices(basis_);
  Eigen::MatrixXd mu(size(), 3);
  for (std::size_t k = 0; k < 3; ++k) {
    const Eigen::MatrixXd mo = C_occ_.transpose() * dipole[k] * C_virt_;
    for (Eigen::Index i = 0; i < n_occ_; ++i)
      for (Eigen::Index a = 0; a < n_virt_; ++a) mu(i * n_virt_ + a, static_cast<Eigen::Index>(k)) = mo(i, a);
  }
  return mu;
}

}  // namespace ellidft
