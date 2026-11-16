#pragma once

#include <array>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "core/molecule.hpp"

namespace ellidft {

// One-electron integral matrices over the basis functions of `basis`
// (McMurchie-Davidson scheme). Every shell pair (a, b) is computed
// separately, including both orders; symmetry of the matrices is not imposed.

// S_mu,nu = <mu|nu>
Eigen::MatrixXd overlap_matrix(const BasisSet& basis);

// T_mu,nu = <mu| -1/2 nabla^2 |nu>
Eigen::MatrixXd kinetic_matrix(const BasisSet& basis);

// V_mu,nu = - sum_A Z_A <mu| 1/|r - R_A| |nu>, point nuclei of `molecule`.
Eigen::MatrixXd nuclear_attraction_matrix(const BasisSet& basis, const Molecule& molecule);

// <mu| (r - origin)_k |nu> for k = x, y, z. No charge factor.
std::array<Eigen::MatrixXd, 3> dipole_matrices(const BasisSet& basis, const Vec3& origin = {0.0, 0.0, 0.0});

// Overlap computed through the Hermite expansion of the shell pairs
// (ShellPair). Same result as overlap_matrix(); used to test the shell pair
// data that the Coulomb integrals are built on.
Eigen::MatrixXd overlap_matrix_from_shell_pairs(const BasisSet& basis);

}  // namespace ellidft
