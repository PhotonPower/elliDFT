#pragma once

#include <vector>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "integrals/shell_pair.hpp"

namespace ellidft {

// Coulomb and exchange matrices for given density matrices:
//
//   J_mu,nu = sum_{lambda sigma} (mu nu | lambda sigma) D_lambda,sigma
//   K_mu,nu = sum_{lambda sigma} (mu lambda | nu sigma) D_lambda,sigma
//
// The density matrices may be arbitrary real matrices, also unsymmetric
// (transition densities in response theory). Implementations are
// interchangeable (direct four-centre integrals, density fitting, GPU);
// SCF and response theory use only this interface (specification, section 4).
class JKBuilder {
 public:
  virtual ~JKBuilder() = default;

  // J and K receive one matrix per density. Either pointer may be null if
  // the matrix is not needed.
  virtual void build(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
                     std::vector<Eigen::MatrixXd>* K) = 0;

  // Long-range exchange for range-separated functionals: K as above with the
  // operator erf(omega r12) / r12 in place of 1 / r12.
  virtual void build_long_range_exchange(const std::vector<Eigen::MatrixXd>& densities, double omega,
                                         std::vector<Eigen::MatrixXd>& K_long_range) = 0;

  // True if the cost of build() drops for small density matrices (screening
  // with the density). The SCF then passes density differences.
  virtual bool screens_with_density() const { return false; }
};

struct DirectJKOptions {
  // Schwarz screening: a shell quartet (ab|cd) is skipped if
  //   sqrt(max (ab|ab)) * sqrt(max (cd|cd)) * max |D|  <  schwarz_threshold.
  // The product of the two square roots bounds every integral of the quartet
  // (Cauchy-Schwarz inequality for the Coulomb inner product). 0 switches
  // the screening off.
  double schwarz_threshold = 1e-12;
};

// Integral-direct construction from four-centre integrals: every unique shell
// quartet that survives the screening is computed once per call and
// contracted with the densities. No integrals are stored. Uses all threads;
// the result depends on the number of threads only through the order of
// summation.
class DirectJK final : public JKBuilder {
 public:
  explicit DirectJK(const BasisSet& basis, const DirectJKOptions& options = {});

  void build(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
             std::vector<Eigen::MatrixXd>* K) override;
  void build_long_range_exchange(const std::vector<Eigen::MatrixXd>& densities, double omega,
                                 std::vector<Eigen::MatrixXd>& K_long_range) override;

  bool screens_with_density() const override { return options_.schwarz_threshold > 0.0; }

  // Shell quartets computed and skipped in the last call of build().
  std::size_t quartets_computed() const { return quartets_computed_; }
  std::size_t quartets_skipped() const { return quartets_skipped_; }

 private:
  void accumulate(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
                  std::vector<Eigen::MatrixXd>* K, double omega);

  DirectJKOptions options_;
  int n_functions_;
  int max_l_;
  std::vector<int> first_function_;
  std::vector<ShellPair> pairs_;    // pairs_[a(a+1)/2 + b] for a >= b
  std::vector<double> schwarz_;     // sqrt(max (ab|ab)) per pair
  std::size_t quartets_computed_ = 0, quartets_skipped_ = 0;
};

}  // namespace ellidft
