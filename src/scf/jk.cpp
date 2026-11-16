#include "scf/jk.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/parallel.hpp"
#include "integrals/two_electron.hpp"

namespace ellidft {

DirectJK::DirectJK(const BasisSet& basis, const DirectJKOptions& options)
    : options_(options), n_functions_(basis.n_functions()), max_l_(basis.max_l()) {
  const auto& shells = basis.shells();
  for (std::size_t s = 0; s < shells.size(); ++s) first_function_.push_back(basis.first_function(s));
  pairs_.reserve(shells.size() * (shells.size() + 1) / 2);
  for (std::size_t a = 0; a < shells.size(); ++a)
    for (std::size_t b = 0; b <= a; ++b) pairs_.push_back(make_shell_pair(shells[a], shells[b]));

  // Schwarz factors: sqrt of the largest diagonal integral (ab|ab) of a pair.
  schwarz_.assign(pairs_.size(), 0.0);
  const auto n_pairs = static_cast<long>(pairs_.size());
#pragma omp parallel
  {
    EriEngine engine(max_l_);
#pragma omp for schedule(dynamic, 16)
    for (long p = 0; p < n_pairs; ++p) {
      const ShellPair& pair = pairs_[static_cast<std::size_t>(p)];
      const std::vector<double>& block = engine.compute(pair, pair);
      const auto nf = static_cast<std::size_t>(pair.nf);
      double largest = 0.0;
      for (std::size_t f = 0; f < nf; ++f) largest = std::max(largest, std::abs(block[f * nf + f]));
      schwarz_[static_cast<std::size_t>(p)] = std::sqrt(largest);
    }
  }
}

void DirectJK::build(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
                     std::vector<Eigen::MatrixXd>* K) {
  accumulate(densities, J, K, 0.0);
}

void DirectJK::build_long_range_exchange(const std::vector<Eigen::MatrixXd>& densities, double omega,
                                         std::vector<Eigen::MatrixXd>& K_long_range) {
  if (!(omega > 0.0)) throw std::invalid_argument("DirectJK: long-range exchange needs omega > 0");
  // The Schwarz factors of the full Coulomb operator also bound the
  // attenuated integrals: erf(omega r)/r <= 1/r pointwise and both kernels
  // are positive definite, so (ab|ab)_erf <= (ab|ab).
  accumulate(densities, nullptr, &K_long_range, omega);
}

void DirectJK::accumulate(const std::vector<Eigen::MatrixXd>& densities, std::vector<Eigen::MatrixXd>* J,
                          std::vector<Eigen::MatrixXd>* K, double omega) {
  const std::size_t n_densities = densities.size();
  double density_max = 0.0;
  for (const Eigen::MatrixXd& D : densities) {
    if (D.rows() != n_functions_ || D.cols() != n_functions_)
      throw std::invalid_argument("DirectJK: density matrix has the wrong dimension");
    if (D.size() > 0) density_max = std::max(density_max, D.cwiseAbs().maxCoeff());
  }
  const bool want_J = J != nullptr, want_K = K != nullptr;

  // J depends only on the symmetric part of D because (mu nu|lambda sigma) is
  // symmetric in lambda, sigma. With S = D + D^T the unique quartets give
  //   Jt_ij += (ij|kl) S_kl,  Jt_kl += (ij|kl) S_ij,  J = Jt + Jt^T.
  std::vector<Eigen::MatrixXd> symmetric_density;
  if (want_J)
    for (const Eigen::MatrixXd& D : densities) symmetric_density.push_back(D + D.transpose());

  const int n_threads = max_threads();
  const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(n_functions_, n_functions_);
  // One set of accumulators per thread, added up in thread order afterwards.
  std::vector<std::vector<Eigen::MatrixXd>> Jt(static_cast<std::size_t>(n_threads)),
      Kt(static_cast<std::size_t>(n_threads));
  std::vector<std::size_t> computed(static_cast<std::size_t>(n_threads), 0), skipped(computed.size(), 0);

  const std::size_t n_shells = first_function_.size();
  const double threshold = options_.schwarz_threshold;

#pragma omp parallel num_threads(n_threads)
  {
    const auto thread = static_cast<std::size_t>(thread_index());
    std::vector<Eigen::MatrixXd>& Jl = Jt[thread];
    std::vector<Eigen::MatrixXd>& Kl = Kt[thread];
    if (want_J) Jl.assign(n_densities, zero);
    if (want_K) Kl.assign(n_densities, zero);
    EriEngine engine(max_l_);

    // Unique shell quartets: a >= b, c >= d, pair (ab) >= pair (cd). The eight
    // permutations (ab|cd) = (ba|cd) = (ab|dc) = (ba|dc) = (cd|ab) = ... are
    // applied explicitly; where shells coincide a permutation maps the block
    // onto itself and the contribution is halved to avoid double counting.
#pragma omp for schedule(dynamic, 1)
    for (long ab_index = static_cast<long>(pairs_.size()) - 1; ab_index >= 0; --ab_index) {
      const auto ab = static_cast<std::size_t>(ab_index);
      // a(a+1)/2 + b = ab with b <= a
      auto a = static_cast<std::size_t>((std::sqrt(8.0 * static_cast<double>(ab) + 1.0) - 1.0) / 2.0);
      while (a * (a + 1) / 2 > ab) --a;
      while ((a + 1) * (a + 2) / 2 <= ab) ++a;
      const std::size_t b = ab - a * (a + 1) / 2;
      const ShellPair& bra = pairs_[ab];
      const double q_ab = schwarz_[ab];

      for (std::size_t c = 0; c <= a && c < n_shells; ++c) {
        for (std::size_t d = 0; d <= c; ++d) {
          const std::size_t cd = c * (c + 1) / 2 + d;
          if (cd > ab) continue;
          if (threshold > 0.0 && q_ab * schwarz_[cd] * density_max < threshold) {
            ++skipped[thread];
            continue;
          }
          ++computed[thread];
          const ShellPair& ket = pairs_[cd];

          double scale = 1.0;
          if (a == b) scale *= 0.5;
          if (c == d) scale *= 0.5;
          if (ab == cd) scale *= 0.5;

          const std::vector<double>& block = engine.compute(bra, ket, omega);
          const int ia = first_function_[a], ib = first_function_[b];
          const int ic = first_function_[c], id = first_function_[d];

          std::size_t index = 0;
          for (int fa = 0; fa < bra.nfa; ++fa) {
            const int i = ia + fa;
            for (int fb = 0; fb < bra.nfb; ++fb) {
              const int j = ib + fb;
              for (int fc = 0; fc < ket.nfa; ++fc) {
                const int k = ic + fc;
                for (int fd = 0; fd < ket.nfb; ++fd, ++index) {
                  const int l = id + fd;
                  const double value = scale * block[index];
                  if (value == 0.0) continue;
                  for (std::size_t n = 0; n < n_densities; ++n) {
                    if (want_J) {
                      Jl[n](i, j) += value * symmetric_density[n](k, l);
                      Jl[n](k, l) += value * symmetric_density[n](i, j);
                    }
                    if (want_K) {
                      // K_mu,nu += (mu lambda|nu sigma) D_lambda,sigma for each
                      // of the eight permutations of (ij|kl).
                      const Eigen::MatrixXd& D = densities[n];
                      Eigen::MatrixXd& Kn = Kl[n];
                      Kn(i, k) += value * D(j, l);  // (ij|kl)
                      Kn(j, k) += value * D(i, l);  // (ji|kl)
                      Kn(i, l) += value * D(j, k);  // (ij|lk)
                      Kn(j, l) += value * D(i, k);  // (ji|lk)
                      Kn(k, i) += value * D(l, j);  // (kl|ij)
                      Kn(l, i) += value * D(k, j);  // (lk|ij)
                      Kn(k, j) += value * D(l, i);  // (kl|ji)
                      Kn(l, j) += value * D(k, i);  // (lk|ji)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  quartets_computed_ = quartets_skipped_ = 0;
  for (std::size_t t = 0; t < computed.size(); ++t) {
    quartets_computed_ += computed[t];
    quartets_skipped_ += skipped[t];
  }
  if (want_J) {
    J->assign(n_densities, zero);
    for (std::size_t n = 0; n < n_densities; ++n) {
      Eigen::MatrixXd sum = zero;
      for (const auto& thread_J : Jt)
        if (!thread_J.empty()) sum += thread_J[n];
      (*J)[n] = sum + sum.transpose();
    }
  }
  if (want_K) {
    K->assign(n_densities, zero);
    for (std::size_t n = 0; n < n_densities; ++n)
      for (const auto& thread_K : Kt)
        if (!thread_K.empty()) (*K)[n] += thread_K[n];
  }
}

}  // namespace ellidft
