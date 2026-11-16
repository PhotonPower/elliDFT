#pragma once

#include <vector>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"

namespace ellidft {

// Values of the basis functions, and optionally their gradients, at a block
// of points: value(p, mu) = chi_mu(r_p); dx, dy, dz likewise for the
// derivatives with respect to the electron coordinate.
struct BasisValues {
  Eigen::MatrixXd value;
  Eigen::MatrixXd dx, dy, dz;  // empty unless gradients were requested
};

// Evaluates all basis functions at the given points.
//
// `cutoff` controls the skipping of distant functions: a shell is left at
// zero for the whole block if an upper bound of its radial part
// |R(r)| r^l at the nearest point of the block is below `cutoff`. With
// cutoff = 0 nothing is skipped. The default is far below the accuracy of
// the grids.
void evaluate_basis(const BasisSet& basis, const Vec3* points, std::size_t n_points, bool gradients,
                    BasisValues& out, double cutoff = 1e-16);

// Same, restricted to the functions that are not skipped in this block:
// `out` has one column per entry of `active`, which lists the indices of
// those basis functions in ascending order. Working with these columns only
// is what makes the quadrature scale with the size of the molecule.
void evaluate_basis_active(const BasisSet& basis, const Vec3* points, std::size_t n_points, bool gradients,
                           BasisValues& out, std::vector<int>& active, double cutoff = 1e-16);

}  // namespace ellidft
