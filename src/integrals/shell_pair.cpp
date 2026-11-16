#include "integrals/shell_pair.hpp"

#include <cmath>

#include <Eigen/Dense>

#include "basis/solid_harmonics.hpp"
#include "integrals/hermite.hpp"

namespace ellidft {

ShellPair make_shell_pair(const Shell& a, const Shell& b) {
  ShellPair pair;
  pair.la = a.l;
  pair.lb = b.l;
  pair.L = a.l + b.l;
  pair.nfa = a.n_functions();
  pair.nfb = b.n_functions();
  pair.nf = pair.nfa * pair.nfb;
  pair.nherm = n_hermite(pair.L);

  const auto& powers_a = cartesian_powers(a.l);
  const auto& powers_b = cartesian_powers(b.l);
  const Eigen::MatrixXd& Ta = shell_transformation(a.l, a.spherical);
  const Eigen::MatrixXd& Tb = shell_transformation(b.l, b.spherical);
  const std::vector<HermiteIndex> hermite = hermite_indices(pair.L);

  const auto nca = static_cast<Eigen::Index>(powers_a.size());
  const auto ncb = static_cast<Eigen::Index>(powers_b.size());
  Eigen::MatrixXd raw(nca, ncb);
  Eigen::MatrixXd transformed(pair.nfa, pair.nfb);

  pair.primitives.reserve(a.n_primitives() * b.n_primitives());
  for (std::size_t ka = 0; ka < a.n_primitives(); ++ka) {
    for (std::size_t kb = 0; kb < b.n_primitives(); ++kb) {
      const double ea = a.exponents[ka];
      const double eb = b.exponents[kb];
      const double coefficient = a.norm_coefficients[ka] * b.norm_coefficients[kb];

      PrimitivePair prim;
      prim.p = ea + eb;
      for (std::size_t k = 0; k < 3; ++k) prim.P[k] = (ea * a.center[k] + eb * b.center[k]) / prim.p;
      prim.E.assign(static_cast<std::size_t>(pair.nherm) * static_cast<std::size_t>(pair.nf), 0.0);

      const HermiteE Ex(a.l, b.l, ea, eb, a.center[0], b.center[0]);
      const HermiteE Ey(a.l, b.l, ea, eb, a.center[1], b.center[1]);
      const HermiteE Ez(a.l, b.l, ea, eb, a.center[2], b.center[2]);

      for (std::size_t h = 0; h < hermite.size(); ++h) {
        const auto [t, u, v] = hermite[h];
        // Product of cartesian Gaussians: E_t^{ax bx} E_u^{ay by} E_v^{az bz}
        // (Hernández Vera 2025, Eq. 46).
        for (Eigen::Index ca = 0; ca < nca; ++ca) {
          const auto& pa = powers_a[static_cast<std::size_t>(ca)];
          for (Eigen::Index cb = 0; cb < ncb; ++cb) {
            const auto& pb = powers_b[static_cast<std::size_t>(cb)];
            raw(ca, cb) = Ex(pa.x, pb.x, t) * Ey(pa.y, pb.y, u) * Ez(pa.z, pb.z, v);
          }
        }
        // Cartesian monomials -> basis functions of both shells (Eq. 62).
        transformed.noalias() = coefficient * (Ta * raw * Tb.transpose());
        double* out = &prim.E[h * static_cast<std::size_t>(pair.nf)];
        for (int fa = 0; fa < pair.nfa; ++fa)
          for (int fb = 0; fb < pair.nfb; ++fb) out[fa * pair.nfb + fb] = transformed(fa, fb);
      }
      pair.primitives.push_back(std::move(prim));
    }
  }
  return pair;
}

ShellPair make_single_shell(const Shell& a) {
  // The constant function 1 as an s shell with exponent 0 at the same centre.
  // The transformation of an s shell multiplies by 1/sqrt(4 pi), which the
  // coefficient sqrt(4 pi) cancels.
  Shell unit;
  unit.atom = a.atom;
  unit.center = a.center;
  unit.l = 0;
  unit.spherical = true;
  unit.exponents = {0.0};
  unit.coefficients = {1.0};
  unit.norm_coefficients = {std::sqrt(4.0 * M_PI)};
  return make_shell_pair(a, unit);
}

}  // namespace ellidft
