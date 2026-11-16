#pragma once

#include <deque>

#include <Eigen/Dense>

namespace ellidft {

// Direct inversion in the iterative subspace for the SCF
// (Hernández Vera 2025, Eqs. 65-69): the next Fock matrix is the linear
// combination sum_i c_i F_i of stored Fock matrices, sum_i c_i = 1, that
// minimises the norm of the combined error matrix sum_i c_i e_i.
class Diis {
 public:
  explicit Diis(std::size_t max_vectors = 8);

  // Stores a Fock matrix with its error matrix and returns the extrapolated
  // Fock matrix. The oldest pair is dropped when max_vectors is exceeded.
  Eigen::MatrixXd extrapolate(const Eigen::MatrixXd& fock, const Eigen::MatrixXd& error);

  std::size_t size() const { return focks_.size(); }

 private:
  std::size_t max_vectors_;
  std::deque<Eigen::MatrixXd> focks_, errors_;
};

}  // namespace ellidft
