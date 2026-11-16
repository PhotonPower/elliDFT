#pragma once

#include <memory>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "core/molecule.hpp"
#include "grid/molecular_grid.hpp"
#include "scf/jk.hpp"
#include "scf/rhf.hpp"
#include "xc/xc_evaluator.hpp"

namespace ellidft {

// The orbital rotation Hessians of linear response theory for a closed-shell
// ground state with real orbitals, singlet adapted (specification, section 8):
//
//   A_{ia,jb} = delta_ij delta_ab (eps_a - eps_i) + 2 (ia|jb) + 2 (ia|f_xc|jb) - c_x (ij|ab)
//   B_{ia,jb} =                                     2 (ia|jb) + 2 (ia|f_xc|jb) - c_x (ib|ja)
//
// i, j occupied, a, b virtual orbitals. For range-separated functionals the
// exchange terms carry the same coefficients and operators as in the Fock
// matrix: c_x (..|..) + c_lr (..|..)_erf.
//
// The matrices are never built. (A+B) and (A-B) are applied to trial vectors
// through the J/K interface and the exchange-correlation kernel: with the
// transition density P = C_occ V C_virt^T,
//
//   (A+B) V = de * V + C_occ^T [ 4 J[P] + 4 K_xc[P] - c_x (K[P] + K[P]^T) ] C_virt
//   (A-B) V = de * V + C_occ^T [                    - c_x (K[P] - K[P]^T) ] C_virt
//
// (and the same with c_lr and K_LR), de the orbital energy differences.
//
// Vectors have n_occ * n_virt entries, index i * n_virt + a.
class ResponseOperator {
 public:
  // `scf` must be the converged result for `molecule`, `basis`, `method`; `jk`
  // must belong to `basis`. The references must outlive the operator.
  ResponseOperator(const Molecule& molecule, const BasisSet& basis, const Method& method, const ScfResult& scf,
                   JKBuilder& jk);

  Eigen::Index n_occ() const { return n_occ_; }
  Eigen::Index n_virt() const { return n_virt_; }
  Eigen::Index size() const { return n_occ_ * n_virt_; }

  // eps_a - eps_i for every pair.
  const Eigen::VectorXd& orbital_energy_differences() const { return differences_; }

  // (A+B) V and (A-B) V for the columns of V. Either result may be null.
  void apply(const Eigen::MatrixXd& V, Eigen::MatrixXd* plus, Eigen::MatrixXd* minus) const;

  // Dipole matrix elements between occupied and virtual orbitals,
  // mu^k_ia = <i| r_k |a>, as three columns (k = x, y, z; origin (0,0,0)).
  Eigen::MatrixXd dipole_vectors() const;

  // Number of calls of apply() and of vectors processed (for reports).
  int applications() const { return applications_; }

 private:
  const BasisSet& basis_;
  Functional functional_;
  JKBuilder& jk_;
  Eigen::MatrixXd C_occ_, C_virt_, density_;
  Eigen::VectorXd differences_;
  Eigen::Index n_occ_, n_virt_;
  std::unique_ptr<XcEvaluator> xc_;
  MolecularGrid grid_;
  mutable int applications_ = 0;
};

}  // namespace ellidft
