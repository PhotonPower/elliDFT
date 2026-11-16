#pragma once

#include <vector>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "grid/molecular_grid.hpp"
#include "xc/xc_evaluator.hpp"

namespace ellidft {

struct XcIntegrals {
  double energy = 0.0;        // E_xc = integral exc(rho, sigma) d^3r
  double n_electrons = 0.0;   // integral rho d^3r on the grid
  Eigen::MatrixXd matrix;     // V_mu,nu = d E_xc / d D_mu,nu (symmetric)
};

// Exchange-correlation energy and matrix for a closed-shell density matrix D
// (rho = sum D_mu,nu chi_mu chi_nu) by quadrature on the molecular grid.
// The grid is processed in blocks of `block_size` points; in every block only
// the basis functions that are not negligible there take part (`basis_cutoff`,
// see evaluate_basis_active; 0 keeps all functions). Only the symmetric part
// of D enters.
XcIntegrals integrate_xc(const BasisSet& basis, const MolecularGrid& grid, const XcEvaluator& functional,
                         const Eigen::MatrixXd& density, std::size_t block_size = 256,
                         double basis_cutoff = 1e-13);

// Exchange-correlation kernel applied to first-order densities:
//
//   K_xc[P]_mu,nu = d V_xc[D]_mu,nu / d D  contracted with P
//                 = sum_{lambda sigma} (mu nu| f_xc |lambda sigma) P_lambda,sigma,
//
// with f_xc the second derivative of the exchange-correlation energy with
// respect to the total density, at the ground-state density matrix
// `density`. One matrix per trial density; only the symmetric part of a
// trial density enters. LDA and GGA.
std::vector<Eigen::MatrixXd> integrate_xc_kernel(const BasisSet& basis, const MolecularGrid& grid,
                                                 const XcEvaluator& functional, const Eigen::MatrixXd& density,
                                                 const std::vector<Eigen::MatrixXd>& trial_densities,
                                                 std::size_t block_size = 256, double basis_cutoff = 1e-13);

}  // namespace ellidft
