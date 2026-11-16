#pragma once

#include <vector>

#include "basis/basis_set.hpp"

namespace ellidft {

// Data of one pair of primitives of a shell pair.
struct PrimitivePair {
  double p = 0.0;   // sum of the exponents
  double P[3]{};    // centre of the product Gaussian
  // Expansion of the products of basis functions in Hermite Gaussians,
  //
  //   chi_fa(r) chi_fb(r)  <-  sum_h E[h * nf + (fa * nfb + fb)] Lambda_h(r - P),
  //
  // restricted to this pair of primitives; h runs over hermite_indices(L).
  // Contraction coefficients, normalisation and the transformation to the
  // basis functions of the shells (spherical or cartesian) are included.
  std::vector<double> E;
};

// A pair of shells prepared for integral evaluation. The products of the
// basis functions of the two shells are expanded in Hermite Gaussians
// (Hernández Vera 2025, Eqs. 46 and 62), primitive pair by primitive pair.
struct ShellPair {
  int la = 0, lb = 0;
  int L = 0;        // la + lb
  int nfa = 0, nfb = 0;
  int nf = 0;       // nfa * nfb
  int nherm = 0;    // n_hermite(L)
  std::vector<PrimitivePair> primitives;
};

ShellPair make_shell_pair(const Shell& a, const Shell& b);

// A single shell in the form of a shell pair: its basis functions expanded in
// Hermite Gaussians (the "second function" is the constant 1). With this,
// EriEngine computes three-centre integrals (ab|P) and two-centre integrals
// (P|Q) of density fitting with the same code as four-centre integrals.
ShellPair make_single_shell(const Shell& a);

}  // namespace ellidft
