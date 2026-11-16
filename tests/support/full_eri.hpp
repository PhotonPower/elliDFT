#pragma once

// Full four-index tensor of two-electron integrals. For tests with small
// molecules only; the library never builds this tensor (specification,
// section 4).

#include <vector>

#include "basis/basis_set.hpp"
#include "integrals/shell_pair.hpp"
#include "integrals/two_electron.hpp"

namespace ellidft::testing {

struct FullEri {
  int n = 0;
  std::vector<double> values;  // n^4 entries
  double operator()(int i, int j, int k, int l) const {
    const auto N = static_cast<std::size_t>(n);
    return values[((static_cast<std::size_t>(i) * N + static_cast<std::size_t>(j)) * N +
                   static_cast<std::size_t>(k)) * N + static_cast<std::size_t>(l)];
  }
  double& at(int i, int j, int k, int l) {
    const auto N = static_cast<std::size_t>(n);
    return values[((static_cast<std::size_t>(i) * N + static_cast<std::size_t>(j)) * N +
                   static_cast<std::size_t>(k)) * N + static_cast<std::size_t>(l)];
  }
};

// Every shell quartet (a, b, c, d) is computed on its own, in all orders; no
// permutational symmetry is used. The symmetry of the result is therefore a
// property of the integral code, not of this function. With omega > 0 the
// operator is erf(omega r12)/r12.
inline FullEri full_eri_tensor(const BasisSet& basis, double omega = 0.0) {
  const auto& shells = basis.shells();
  const std::size_t ns = shells.size();
  std::vector<ShellPair> pairs;
  pairs.reserve(ns * ns);
  for (std::size_t a = 0; a < ns; ++a)
    for (std::size_t b = 0; b < ns; ++b) pairs.push_back(make_shell_pair(shells[a], shells[b]));

  FullEri eri;
  eri.n = basis.n_functions();
  const auto N = static_cast<std::size_t>(eri.n);
  eri.values.assign(N * N * N * N, 0.0);

  EriEngine engine(basis.max_l());
  for (std::size_t a = 0; a < ns; ++a)
    for (std::size_t b = 0; b < ns; ++b)
      for (std::size_t c = 0; c < ns; ++c)
        for (std::size_t d = 0; d < ns; ++d) {
          const ShellPair& bra = pairs[a * ns + b];
          const ShellPair& ket = pairs[c * ns + d];
          const std::vector<double>& block = engine.compute(bra, ket, omega);
          const int ia = basis.first_function(a), ib = basis.first_function(b);
          const int ic = basis.first_function(c), id = basis.first_function(d);
          std::size_t k = 0;
          for (int fa = 0; fa < bra.nfa; ++fa)
            for (int fb = 0; fb < bra.nfb; ++fb)
              for (int fc = 0; fc < ket.nfa; ++fc)
                for (int fd = 0; fd < ket.nfb; ++fd) eri.at(ia + fa, ib + fb, ic + fc, id + fd) = block[k++];
        }
  return eri;
}

}  // namespace ellidft::testing
