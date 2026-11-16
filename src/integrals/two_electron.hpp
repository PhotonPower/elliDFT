#pragma once

#include <vector>

#include "integrals/hermite.hpp"
#include "integrals/shell_pair.hpp"

namespace ellidft {

// Two-electron repulsion integrals (ab|cd) over the basis functions of two
// shell pairs, McMurchie-Davidson scheme (Hernández Vera 2025, Eq. 59).
//
// The engine owns its work arrays; use one engine per thread.
class EriEngine {
 public:
  // `max_l`: highest angular momentum of any shell that will be passed in.
  explicit EriEngine(int max_l);

  // Integrals for all basis functions of the four shells. The result has
  // bra.nf * ket.nf entries with index
  //
  //   ((fa * nfb + fb) * nfc + fc) * nfd + fd
  //
  // and stays valid until the next call. Chemists' notation: functions a, b
  // belong to electron 1, c, d to electron 2.
  //
  // With omega > 0 the operator is erf(omega r12) / r12 instead of 1 / r12
  // (long-range part of the Coulomb interaction, in 1/Bohr).
  const std::vector<double>& compute(const ShellPair& bra, const ShellPair& ket, double omega = 0.0);

 private:
  int max_L_pair_;  // 2 * max_l
  HermiteCoulomb R_;
  std::vector<std::vector<HermiteIndex>> hermite_;  // hermite_[L] = hermite_indices(L)
  std::vector<double> result_, half_;
  std::vector<std::size_t> bra_offset_, ket_offset_;
  std::vector<double> ket_sign_;
};

}  // namespace ellidft
