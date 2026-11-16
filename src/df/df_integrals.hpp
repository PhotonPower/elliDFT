#pragma once

#include <Eigen/Dense>

#include "basis/basis_set.hpp"

namespace ellidft {

// Two-centre Coulomb integrals of the auxiliary basis, the metric of density
// fitting: V_PQ = (P|Q) = integral P(r1) Q(r2) / |r1 - r2|.
Eigen::MatrixXd two_center_coulomb(const BasisSet& auxiliary, double omega = 0.0);

// Three-centre Coulomb integrals (mu nu|P) for mu >= nu: row index
// mu (mu + 1) / 2 + nu, column index P.
Eigen::MatrixXd three_center_coulomb(const BasisSet& basis, const BasisSet& auxiliary, double omega = 0.0);

}  // namespace ellidft
