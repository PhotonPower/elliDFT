#include <cmath>
#include <stdexcept>

#include "xc/xc_evaluator.hpp"

namespace ellidft {

namespace {

// Densities below this value give zero energy and potential.
constexpr double density_threshold = 1e-24;

// VWN5 correlation, paramagnetic (spin-unpolarised) case. Parameters and
// formula as in libxc 7.0.0, maple/vwn.mpl (A_vwn[1], b_vwn[1], c_vwn[1],
// x0_vwn[1] and f_aux), consulted 2026-11-14; energies in Hartree.
constexpr double vwn_A = 0.0310907;
constexpr double vwn_b = 3.72744;
constexpr double vwn_c = 12.9352;
constexpr double vwn_x0 = -0.10498;

class OwnLda final : public XcEvaluator {
 public:
  bool needs_gradient() const override { return false; }

  void evaluate(std::size_t n, const double* rho, const double* /*sigma*/, double* exc, double* vrho,
                double* vsigma) const override {
    const double Q = std::sqrt(4.0 * vwn_c - vwn_b * vwn_b);
    const double f1 = 2.0 * vwn_b / Q;
    const double f2 = vwn_b * vwn_x0 / (vwn_x0 * vwn_x0 + vwn_b * vwn_x0 + vwn_c);
    const double f3 = 2.0 * (2.0 * vwn_x0 + vwn_b) / Q;
    // Slater exchange per particle: -3/4 (3/pi)^(1/3) rho^(1/3)
    // (Wikipedia, "Local-density approximation", exchange functional in
    // atomic units, consulted 2026-11-14).
    const double cx = -0.75 * std::cbrt(3.0 / M_PI);

    for (std::size_t i = 0; i < n; ++i) {
      if (vsigma != nullptr) vsigma[i] = 0.0;
      const double r = rho[i];
      if (!(r > density_threshold)) {
        exc[i] = 0.0;
        vrho[i] = 0.0;
        continue;
      }
      // Exchange: e_x = cx rho^(1/3), v_x = d(rho e_x)/d rho = 4/3 e_x.
      const double ex = cx * std::cbrt(r);
      const double vx = 4.0 / 3.0 * ex;

      // Correlation. Wigner-Seitz radius: 4 pi rs^3 / 3 = 1 / rho.
      const double rs = std::cbrt(3.0 / (4.0 * M_PI * r));
      const double x = std::sqrt(rs);
      const double X = rs + vwn_b * x + vwn_c;
      const double atan_term = std::atan(Q / (2.0 * x + vwn_b));
      // f_aux of vwn.mpl:
      //   e_c = A [ ln(rs/X) + (f1 - f2 f3) atan(Q/(2x+b)) - f2 ln((x-x0)^2/X) ]
      const double ec = vwn_A * (std::log(rs / X) + (f1 - f2 * f3) * atan_term -
                                 f2 * std::log((x - vwn_x0) * (x - vwn_x0) / X));
      // Derivative with respect to rs, term by term, with dX/drs = 1 + b/(2x)
      // and (2x+b)^2 + Q^2 = 4X:
      //   d/drs ln(rs/X)            = 1/rs - X'/X
      //   d/drs atan(Q/(2x+b))      = -Q / (4 x X)
      //   d/drs ln((x-x0)^2/X)      = 1/(x (x-x0)) - X'/X
      const double dX = 1.0 + vwn_b / (2.0 * x);
      const double dec_drs = vwn_A * (1.0 / rs - dX / X - (f1 - f2 * f3) * Q / (4.0 * x * X) -
                                      f2 * (1.0 / (x * (x - vwn_x0)) - dX / X));
      // d rs / d rho = -rs / (3 rho), hence v_c = e_c - (rs/3) de_c/drs.
      const double vc = ec - rs / 3.0 * dec_drs;

      exc[i] = r * (ex + ec);
      vrho[i] = vx + vc;
    }
  }

  void evaluate_kernel(std::size_t, const double*, const double*, double*, double*, double*, double*,
                       double*) const override {
    throw std::logic_error("the own LDA does not provide second derivatives; use the libxc evaluator");
  }
};

}  // namespace

std::unique_ptr<XcEvaluator> make_own_lda_evaluator() { return std::make_unique<OwnLda>(); }

}  // namespace ellidft
