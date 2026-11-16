#include "integrals/two_electron.hpp"

#include <cmath>
#include <stdexcept>

namespace ellidft {

EriEngine::EriEngine(int max_l) : max_L_pair_(2 * max_l), R_(4 * max_l) {
  for (int L = 0; L <= max_L_pair_; ++L) hermite_.push_back(hermite_indices(L));
}

const std::vector<double>& EriEngine::compute(const ShellPair& bra, const ShellPair& ket, double omega) {
  if (bra.L > max_L_pair_ || ket.L > max_L_pair_)
    throw std::invalid_argument("EriEngine: shell pair exceeds the angular momentum of the engine");

  const int L = bra.L + ket.L;
  const auto dim = static_cast<std::size_t>(L + 1);
  const auto nhb = static_cast<std::size_t>(bra.nherm);
  const auto nhk = static_cast<std::size_t>(ket.nherm);
  const auto nfab = static_cast<std::size_t>(bra.nf);
  const auto nfcd = static_cast<std::size_t>(ket.nf);

  // Position of R_{tuv} in the cube of HermiteCoulomb is (t*dim + u)*dim + v,
  // so the position of R_{t+tau,u+nu,v+phi} is the sum of two offsets.
  const auto& hb_list = hermite_[static_cast<std::size_t>(bra.L)];
  const auto& hk_list = hermite_[static_cast<std::size_t>(ket.L)];
  bra_offset_.resize(nhb);
  ket_offset_.resize(nhk);
  ket_sign_.resize(nhk);
  for (std::size_t h = 0; h < nhb; ++h)
    bra_offset_[h] = (static_cast<std::size_t>(hb_list[h].t) * dim + static_cast<std::size_t>(hb_list[h].u)) * dim +
                     static_cast<std::size_t>(hb_list[h].v);
  for (std::size_t h = 0; h < nhk; ++h) {
    ket_offset_[h] = (static_cast<std::size_t>(hk_list[h].t) * dim + static_cast<std::size_t>(hk_list[h].u)) * dim +
                     static_cast<std::size_t>(hk_list[h].v);
    // (-1)^(tau + nu + phi), Eq. 59
    ket_sign_[h] = ((hk_list[h].t + hk_list[h].u + hk_list[h].v) % 2 == 0) ? 1.0 : -1.0;
  }

  result_.assign(nfab * nfcd, 0.0);
  half_.resize(nhb * nfcd);
  const double two_pi_52 = 2.0 * std::pow(M_PI, 2.5);

  // Hernández Vera 2025, Eq. 59, for every pair of primitive pairs:
  //   [ab|cd] = 2 pi^(5/2) / (p q sqrt(p+q)) sum_{tuv} sum_{tau nu phi} (-1)^(tau+nu+phi)
  //             E^{ab}_{tuv} E^{cd}_{tau nu phi} R^0_{t+tau,u+nu,v+phi}(alpha, P - Q),
  // alpha = p q / (p + q). The sum over the ket primitives and ket Hermite
  // indices is done first (intermediate of Eq. 64), then the bra is contracted.
  for (const PrimitivePair& bp : bra.primitives) {
    std::fill(half_.begin(), half_.end(), 0.0);
    for (const PrimitivePair& kp : ket.primitives) {
      const double p = bp.p, q = kp.p;
      double alpha = p * q / (p + q);
      double prefactor = two_pi_52 / (p * q * std::sqrt(p + q));
      if (omega > 0.0) {
        // Operator erf(omega r)/r: only the Boys function changes
        // (specification, section 5),
        //   F_n(T) -> s^(n+1/2) F_n(s T),  s = omega^2 / (omega^2 + alpha).
        // The source terms (-2 alpha)^n F_n(alpha R^2) of the Hermite Coulomb
        // integrals therefore become sqrt(s) (-2 s alpha)^n F_n(s alpha R^2):
        // the reduced exponent is scaled by s and the prefactor by sqrt(s).
        const double s = omega * omega / (omega * omega + alpha);
        alpha *= s;
        prefactor *= std::sqrt(s);
      }
      R_.compute(L, alpha, bp.P[0] - kp.P[0], bp.P[1] - kp.P[1], bp.P[2] - kp.P[2]);
      const double* r = R_.data();
      for (std::size_t hb = 0; hb < nhb; ++hb) {
        double* half_row = &half_[hb * nfcd];
        const double* r_row = r + bra_offset_[hb];
        for (std::size_t hk = 0; hk < nhk; ++hk) {
          const double w = prefactor * ket_sign_[hk] * r_row[ket_offset_[hk]];
          const double* e = &kp.E[hk * nfcd];
          for (std::size_t f = 0; f < nfcd; ++f) half_row[f] += w * e[f];
        }
      }
    }
    for (std::size_t hb = 0; hb < nhb; ++hb) {
      const double* half_row = &half_[hb * nfcd];
      const double* e = &bp.E[hb * nfab];
      for (std::size_t fab = 0; fab < nfab; ++fab) {
        const double x = e[fab];
        if (x == 0.0) continue;
        double* out = &result_[fab * nfcd];
        for (std::size_t f = 0; f < nfcd; ++f) out[f] += x * half_row[f];
      }
    }
  }
  return result_;
}

}  // namespace ellidft
