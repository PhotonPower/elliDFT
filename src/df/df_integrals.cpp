#include "df/df_integrals.hpp"

#include <algorithm>
#include <vector>

#include "integrals/shell_pair.hpp"
#include "integrals/two_electron.hpp"

namespace ellidft {

Eigen::MatrixXd two_center_coulomb(const BasisSet& auxiliary, double omega) {
  const auto& shells = auxiliary.shells();
  const Eigen::Index n = auxiliary.n_functions();
  Eigen::MatrixXd V(n, n);
  std::vector<ShellPair> singles;
  singles.reserve(shells.size());
  for (const Shell& s : shells) singles.push_back(make_single_shell(s));

  const auto n_shells = static_cast<long>(shells.size());
#pragma omp parallel
  {
    EriEngine engine(auxiliary.max_l());
#pragma omp for schedule(dynamic, 4)
    for (long p = 0; p < n_shells; ++p) {
      const ShellPair& bra = singles[static_cast<std::size_t>(p)];
      const int ip = auxiliary.first_function(static_cast<std::size_t>(p));
      for (long q = 0; q <= p; ++q) {
        const ShellPair& ket = singles[static_cast<std::size_t>(q)];
        const int iq = auxiliary.first_function(static_cast<std::size_t>(q));
        const std::vector<double>& block = engine.compute(bra, ket, omega);
        for (int fp = 0; fp < bra.nf; ++fp)
          for (int fq = 0; fq < ket.nf; ++fq) {
            const double value = block[static_cast<std::size_t>(fp * ket.nf + fq)];
            V(ip + fp, iq + fq) = value;
            V(iq + fq, ip + fp) = value;
          }
      }
    }
  }
  return V;
}

Eigen::MatrixXd three_center_coulomb(const BasisSet& basis, const BasisSet& auxiliary, double omega) {
  const auto& shells = basis.shells();
  const auto& aux_shells = auxiliary.shells();
  const Eigen::Index n = basis.n_functions();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n * (n + 1) / 2, auxiliary.n_functions());

  std::vector<ShellPair> singles;
  singles.reserve(aux_shells.size());
  for (const Shell& s : aux_shells) singles.push_back(make_single_shell(s));

  const int max_l = std::max(basis.max_l(), auxiliary.max_l());
  const auto n_shells = static_cast<long>(shells.size());
#pragma omp parallel
  {
    EriEngine engine(max_l);
#pragma omp for schedule(dynamic, 1)
    for (long a = n_shells - 1; a >= 0; --a) {
      for (long b = 0; b <= a; ++b) {
        const ShellPair bra = make_shell_pair(shells[static_cast<std::size_t>(a)], shells[static_cast<std::size_t>(b)]);
        const int ia = basis.first_function(static_cast<std::size_t>(a));
        const int ib = basis.first_function(static_cast<std::size_t>(b));
        for (std::size_t p = 0; p < aux_shells.size(); ++p) {
          const ShellPair& ket = singles[p];
          const int ip = auxiliary.first_function(p);
          const std::vector<double>& block = engine.compute(bra, ket, omega);
          for (int fa = 0; fa < bra.nfa; ++fa)
            for (int fb = 0; fb < bra.nfb; ++fb) {
              const int i = ia + fa, j = ib + fb;
              if (j > i) continue;  // same shell: keep the lower triangle
              const Eigen::Index row = static_cast<Eigen::Index>(i) * (i + 1) / 2 + j;
              const double* values = &block[static_cast<std::size_t>((fa * bra.nfb + fb) * ket.nf)];
              for (int fp = 0; fp < ket.nf; ++fp) A(row, ip + fp) = values[fp];
            }
        }
      }
    }
  }
  return A;
}

}  // namespace ellidft
