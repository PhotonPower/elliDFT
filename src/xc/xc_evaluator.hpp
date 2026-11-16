#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "xc/functional.hpp"

namespace ellidft {

// Evaluates the density functional part of a method at grid points, for a
// closed-shell (spin-unpolarised) density.
//
// Inputs per point: the total density rho and, for gradient-corrected
// functionals, sigma = |grad rho|^2.
// Outputs per point:
//   exc     energy per unit volume, E_xc = integral exc d^3r
//   vrho    d exc / d rho    at fixed sigma
//   vsigma  d exc / d sigma  at fixed rho   (zero for LDA)
class XcEvaluator {
 public:
  virtual ~XcEvaluator() = default;
  virtual bool needs_gradient() const = 0;

  // `sigma` may be null if needs_gradient() is false; `vsigma` may be null
  // in that case too, otherwise it is filled (with zeros for LDA).
  virtual void evaluate(std::size_t n, const double* rho, const double* sigma, double* exc, double* vrho,
                        double* vsigma) const = 0;

  // First and second derivatives of the energy per volume, for the
  // exchange-correlation kernel of response theory:
  //   vrho = d exc/d rho, vsigma = d exc/d sigma,
  //   v2rho2 = d2 exc/d rho2, v2rhosigma = d2 exc/(d rho d sigma),
  //   v2sigma2 = d2 exc/d sigma2.
  // For LDA the sigma derivatives are zero. Implementations without second
  // derivatives throw std::logic_error.
  virtual void evaluate_kernel(std::size_t n, const double* rho, const double* sigma, double* vrho, double* vsigma,
                               double* v2rho2, double* v2rhosigma, double* v2sigma2) const = 0;
};

// All components of a Functional through libxc. This is the only place that
// touches libxc (specification, section 2: libxc behind an own interface).
std::unique_ptr<XcEvaluator> make_libxc_evaluator(const Functional& functional);

// Version string of the linked libxc.
const char* libxc_version();

// Exact-exchange parameters of a hybrid functional of libxc, converted to the
// form used in Functional:
//
//   exchange operator = exact_exchange * K + long_range_exchange * K_LR(omega).
//
// libxc reports (omega, alpha, beta) with alpha the coefficient of full exact
// exchange and beta the coefficient of short-range exact exchange,
// alpha K + beta K_SR, K_SR = K - K_LR. Hence exact_exchange = alpha + beta
// and long_range_exchange = -beta. For CAM-B3LYP this gives 0.19 and 0.46
// with omega = 0.33. Functionals without exact exchange give zeros.
struct HybridParameters {
  double exact_exchange = 0.0;
  double long_range_exchange = 0.0;
  double omega = 0.0;
  bool needs_gradient = false;
};
HybridParameters libxc_hybrid_parameters(const std::string& libxc_name);

// Own implementation of the LDA: Slater exchange plus VWN5 correlation, as a
// cross-check of the libxc path (specification, section 2).
std::unique_ptr<XcEvaluator> make_own_lda_evaluator();

}  // namespace ellidft
