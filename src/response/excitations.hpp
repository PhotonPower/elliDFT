#pragma once

#include <vector>

#include <Eigen/Dense>

#include "core/molecule.hpp"
#include "response/response_operator.hpp"

namespace ellidft {

struct ExcitationOptions {
  int n_states = 10;
  // false: full linear response (RPA / Casida); true: Tamm-Dancoff
  // approximation (B = 0). Only singlet excitations are computed.
  bool tamm_dancoff = false;
  double tolerance = 1e-6;   // norm of the residual vectors
  int max_iterations = 100;
};

struct Excitation {
  double energy = 0.0;               // Hartree
  double oscillator_strength = 0.0;  // f = 2/3 omega |t|^2
  Vec3 transition_dipole{};          // length representation, e * Bohr
  Eigen::VectorXd x_plus_y;          // X + Y, normalised to (X+Y)^T (X-Y) = 1
  Eigen::VectorXd x_minus_y;         // X - Y
  double residual = 0.0;
};

struct ExcitationResult {
  bool tamm_dancoff = false;
  bool converged = false;
  int iterations = 0;
  std::vector<Excitation> states;    // ascending energy
};

// Lowest singlet excitations (specification, section 8):
//
//   (A B; B A) (X; Y) = omega (1 0; 0 -1) (X; Y),   (X+Y)^T (X-Y) = 1,
//   t^k = sqrt(2) sum_ia mu^k_ia (X+Y)_ia,   f = 2/3 omega sum_k (t^k)^2.
//
// In terms of Z = X + Y and W = X - Y: (A+B) Z = omega W, (A-B) W = omega Z.
// Solved by a Davidson method in a growing subspace. In the Tamm-Dancoff
// approximation A X = omega X with X^T X = 1 and t = sqrt(2) mu.X.
//
// Throws std::runtime_error if the ground state is unstable (A-B or the
// reduced problem not positive definite).
ExcitationResult excitations(const ResponseOperator& op, const ExcitationOptions& options = {});

}  // namespace ellidft
