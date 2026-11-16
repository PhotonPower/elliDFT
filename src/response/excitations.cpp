#include "response/excitations.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ellidft {

namespace {

// Appends candidate vectors to an orthonormal basis; see polarizability.cpp.
Eigen::Index extend_basis(Eigen::MatrixXd& basis, const Eigen::MatrixXd& candidates, double threshold = 1e-8) {
  Eigen::Index added = 0;
  for (Eigen::Index c = 0; c < candidates.cols(); ++c) {
    Eigen::VectorXd v = candidates.col(c);
    const double initial = v.norm();
    if (!(initial > 0.0)) continue;
    v /= initial;
    for (int pass = 0; pass < 2; ++pass)
      if (basis.cols() > 0) v -= basis * (basis.transpose() * v);
    const double norm = v.norm();
    if (norm < threshold) continue;
    basis.conservativeResize(candidates.rows(), basis.cols() + 1);
    basis.col(basis.cols() - 1) = v / norm;
    ++added;
  }
  return added;
}

Eigen::MatrixXd symmetrized(const Eigen::MatrixXd& m) { return 0.5 * (m + m.transpose()); }

}  // namespace

ExcitationResult excitations(const ResponseOperator& op, const ExcitationOptions& options) {
  const Eigen::Index N = op.size();
  const Eigen::Index n_states = std::min<Eigen::Index>(options.n_states, N);
  if (n_states <= 0) throw std::invalid_argument("excitations: at least one state must be requested");
  const Eigen::VectorXd& d = op.orbital_energy_differences();
  const Eigen::MatrixXd mu = op.dipole_vectors();
  const bool tda = options.tamm_dancoff;

  ExcitationResult result;
  result.tamm_dancoff = tda;

  // Start vectors: unit vectors of the orbital pairs with the smallest energy
  // differences.
  std::vector<Eigen::Index> order(static_cast<std::size_t>(N));
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](Eigen::Index a, Eigen::Index b) { return d(a) < d(b); });
  const Eigen::Index n_start = std::min<Eigen::Index>(N, n_states + std::max<Eigen::Index>(4, n_states / 2));
  Eigen::MatrixXd basis = Eigen::MatrixXd::Zero(N, n_start);
  for (Eigen::Index k = 0; k < n_start; ++k) basis(order[static_cast<std::size_t>(k)], k) = 1.0;

  Eigen::MatrixXd plus(N, 0), minus(N, 0);  // (A+B) b and (A-B) b
  Eigen::Index processed = 0;
  Eigen::VectorXd omega(n_states), residual_norm(n_states);
  Eigen::MatrixXd Z(N, n_states), W(N, n_states);

  for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
    result.iterations = iteration;
    Eigen::MatrixXd new_plus, new_minus;
    op.apply(basis.rightCols(basis.cols() - processed), &new_plus, &new_minus);
    plus.conservativeResize(N, basis.cols());
    minus.conservativeResize(N, basis.cols());
    plus.rightCols(new_plus.cols()) = new_plus;
    minus.rightCols(new_minus.cols()) = new_minus;
    processed = basis.cols();

    Eigen::MatrixXd corrections(N, 0);
    auto add_correction = [&](const Eigen::VectorXd& t) {
      corrections.conservativeResize(N, corrections.cols() + 1);
      corrections.col(corrections.cols() - 1) = t;
    };

    if (tda) {
      // A = [(A+B) + (A-B)] / 2. Reduced problem: (b^T A b) zeta = omega zeta.
      const Eigen::MatrixXd Ab = 0.5 * (plus + minus);
      const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetrized(basis.transpose() * Ab));
      for (Eigen::Index s = 0; s < n_states; ++s) {
        const Eigen::VectorXd zeta = eigen.eigenvectors().col(s);
        omega(s) = eigen.eigenvalues()(s);
        Z.col(s) = basis * zeta;  // X, normalised to 1
        W.col(s) = Z.col(s);
        const Eigen::VectorXd r = Ab * zeta - omega(s) * Z.col(s);
        residual_norm(s) = r.norm();
        if (residual_norm(s) < options.tolerance) continue;
        Eigen::VectorXd t(N);
        for (Eigen::Index k = 0; k < N; ++k) {
          double denominator = d(k) - omega(s);
          if (std::abs(denominator) < 1e-6) denominator = 1e-6;
          t(k) = r(k) / denominator;
        }
        add_correction(t);
      }
    } else {
      // Reduced problem with Sp = b^T (A+B) b and Sm = b^T (A-B) b:
      //   Sm Sp zeta = omega^2 zeta. With zeta = Sm^(1/2) y this is the
      // symmetric problem Sm^(1/2) Sp Sm^(1/2) y = omega^2 y.
      const Eigen::MatrixXd Sp = symmetrized(basis.transpose() * plus);
      const Eigen::MatrixXd Sm = symmetrized(basis.transpose() * minus);
      const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> minus_eigen(Sm);
      if (minus_eigen.eigenvalues().minCoeff() <= 0.0)
        throw std::runtime_error("excitations: A - B is not positive definite (unstable ground state)");
      const Eigen::MatrixXd Sm_root = minus_eigen.eigenvectors() *
                                      minus_eigen.eigenvalues().cwiseSqrt().asDiagonal() *
                                      minus_eigen.eigenvectors().transpose();
      const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(symmetrized(Sm_root * Sp * Sm_root));
      for (Eigen::Index s = 0; s < n_states; ++s) {
        const double omega2 = eigen.eigenvalues()(s);
        if (omega2 <= 0.0)
          throw std::runtime_error("excitations: negative squared excitation energy (unstable ground state)");
        omega(s) = std::sqrt(omega2);
        // Normalisation (X+Y)^T (X-Y) = zeta^T Sp zeta / omega = 1.
        const Eigen::VectorXd zeta = Sm_root * eigen.eigenvectors().col(s) / std::sqrt(omega(s));
        const Eigen::VectorXd w = Sp * zeta / omega(s);
        Z.col(s) = basis * zeta;
        W.col(s) = basis * w;
        // Residuals of (A+B) Z = omega W and (A-B) W = omega Z.
        const Eigen::VectorXd r1 = plus * zeta - omega(s) * W.col(s);
        const Eigen::VectorXd r2 = minus * w - omega(s) * Z.col(s);
        residual_norm(s) = std::max(r1.norm(), r2.norm());
        if (residual_norm(s) < options.tolerance) continue;
        // Preconditioner: A+B and A-B replaced by the orbital energy
        // differences, a 2x2 system per orbital pair.
        Eigen::VectorXd t1(N), t2(N);
        for (Eigen::Index k = 0; k < N; ++k) {
          double denominator = d(k) * d(k) - omega2;
          if (std::abs(denominator) < 1e-6) denominator = 1e-6;
          t1(k) = (d(k) * r1(k) + omega(s) * r2(k)) / denominator;
          t2(k) = (d(k) * r2(k) + omega(s) * r1(k)) / denominator;
        }
        add_correction(t1);
        add_correction(t2);
      }
    }

    if (corrections.cols() == 0) {
      result.converged = true;
      break;
    }
    if (extend_basis(basis, corrections) == 0 || basis.cols() > N) break;
  }

  for (Eigen::Index s = 0; s < n_states; ++s) {
    Excitation state;
    state.energy = omega(s);
    state.x_plus_y = Z.col(s);
    state.x_minus_y = W.col(s);
    state.residual = residual_norm(s);
    double strength = 0.0;
    for (Eigen::Index k = 0; k < 3; ++k) {
      // t^k = sqrt(2) sum_ia mu^k_ia (X+Y)_ia
      const double t = std::sqrt(2.0) * mu.col(k).dot(Z.col(s));
      state.transition_dipole[static_cast<std::size_t>(k)] = t;
      strength += t * t;
    }
    state.oscillator_strength = 2.0 / 3.0 * state.energy * strength;
    result.states.push_back(std::move(state));
  }
  return result;
}

}  // namespace ellidft
