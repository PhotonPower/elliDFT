#pragma once

#include <Eigen/Dense>

#include "response/response_operator.hpp"

namespace ellidft {

struct LinearResponseOptions {
  // Converged when the residual of every equation, relative to the norm of
  // its right-hand side, is below this value.
  double tolerance = 1e-8;
  int max_iterations = 100;
};

struct Polarizability {
  double omega = 0.0;     // frequency, Hartree
  double gamma = 0.0;     // damping, Hartree
  Eigen::Matrix3d real = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d imag = Eigen::Matrix3d::Zero();  // zero without damping
  bool converged = false;
  int iterations = 0;
  double residual = 0.0;  // largest relative residual at the end
};

// Polarisability tensor alpha(omega + i gamma) in atomic units from the
// linear response equations (specification, section 8):
//
//   [ (A B; B A) - (omega + i gamma) (1 0; 0 -1) ] (X; Y) = (mu; mu),
//   alpha_kl = 2 sum_ia mu^k_ia (X^l + Y^l)_ia.
//
// In terms of Z = X + Y and W = X - Y, with z = omega + i gamma:
//
//   (A+B) Z - z W = 2 mu,   (A-B) W - z Z = 0.
//
// For omega = gamma = 0 this is the coupled-perturbed (CPKS) equation
// (A+B) Z = 2 mu. The equations are solved in a growing subspace of trial
// vectors; only products of (A+B) and (A-B) with vectors are needed.
Polarizability polarizability(const ResponseOperator& op, double omega = 0.0, double gamma = 0.0,
                              const LinearResponseOptions& options = {});

}  // namespace ellidft
