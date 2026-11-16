#include <stdexcept>
#include <vector>

#include <xc.h>

#include "xc/xc_evaluator.hpp"

namespace ellidft {

namespace {

// libxc conventions for the unpolarised case (checked by the unit tests
// through numerical differentiation): zk is the energy per particle, so the
// energy per volume is rho * zk; vrho = d(rho zk)/d rho and
// vsigma = d(rho zk)/d sigma with sigma = |grad rho|^2 of the total density.
class LibxcEvaluator final : public XcEvaluator {
 public:
  explicit LibxcEvaluator(const Functional& functional) : needs_gradient_(false) {
    for (const Functional::Component& component : functional.components) {
      const int id = xc_functional_get_number(component.libxc_name.c_str());
      if (id <= 0) {
        release();
        throw std::invalid_argument("libxc does not know the functional '" + component.libxc_name + "'");
      }
      Term term;
      term.coefficient = component.coefficient;
      if (xc_func_init(&term.func, id, XC_UNPOLARIZED) != 0) {
        release();
        throw std::runtime_error("libxc could not initialise '" + component.libxc_name + "'");
      }
      const int family = term.func.info->family;
      if (family == XC_FAMILY_LDA) {
        term.gga = false;
      } else if (family == XC_FAMILY_GGA || family == XC_FAMILY_HYB_GGA) {
        // For hybrids libxc evaluates the density functional part only; the
        // exact exchange is added by the SCF (see libxc_hybrid_parameters).
        term.gga = true;
        needs_gradient_ = true;
      } else {
        xc_func_end(&term.func);
        release();
        throw std::invalid_argument("functional '" + component.libxc_name +
                                    "' is neither LDA nor GGA; other families are not supported");
      }
      terms_.push_back(term);
    }
  }

  ~LibxcEvaluator() override { release(); }
  LibxcEvaluator(const LibxcEvaluator&) = delete;
  LibxcEvaluator& operator=(const LibxcEvaluator&) = delete;

  bool needs_gradient() const override { return needs_gradient_; }

  void evaluate(std::size_t n, const double* rho, const double* sigma, double* exc, double* vrho,
                double* vsigma) const override {
    if (needs_gradient_ && (sigma == nullptr || vsigma == nullptr))
      throw std::invalid_argument("gradient-corrected functional needs sigma and vsigma");
    for (std::size_t i = 0; i < n; ++i) {
      exc[i] = 0.0;
      vrho[i] = 0.0;
      if (vsigma != nullptr) vsigma[i] = 0.0;
    }
    std::vector<double> zk(n), v1(n), v2(n);
    for (const Term& term : terms_) {
      if (term.gga) {
        xc_gga_exc_vxc(&term.func, n, rho, sigma, zk.data(), v1.data(), v2.data());
        for (std::size_t i = 0; i < n; ++i) vsigma[i] += term.coefficient * v2[i];
      } else {
        xc_lda_exc_vxc(&term.func, n, rho, zk.data(), v1.data());
      }
      for (std::size_t i = 0; i < n; ++i) {
        exc[i] += term.coefficient * rho[i] * zk[i];
        vrho[i] += term.coefficient * v1[i];
      }
    }
  }

  void evaluate_kernel(std::size_t n, const double* rho, const double* sigma, double* vrho, double* vsigma,
                       double* v2rho2, double* v2rhosigma, double* v2sigma2) const override {
    if (needs_gradient_ && sigma == nullptr)
      throw std::invalid_argument("gradient-corrected functional needs sigma");
    for (std::size_t i = 0; i < n; ++i) vrho[i] = vsigma[i] = v2rho2[i] = v2rhosigma[i] = v2sigma2[i] = 0.0;
    std::vector<double> zk(n), v1(n), v2(n), f1(n), f2(n), f3(n);
    for (const Term& term : terms_) {
      if (term.gga) {
        xc_gga_exc_vxc_fxc(&term.func, n, rho, sigma, zk.data(), v1.data(), v2.data(), f1.data(), f2.data(),
                           f3.data());
        for (std::size_t i = 0; i < n; ++i) {
          vsigma[i] += term.coefficient * v2[i];
          v2rhosigma[i] += term.coefficient * f2[i];
          v2sigma2[i] += term.coefficient * f3[i];
        }
      } else {
        xc_lda_exc_vxc_fxc(&term.func, n, rho, zk.data(), v1.data(), f1.data());
      }
      for (std::size_t i = 0; i < n; ++i) {
        vrho[i] += term.coefficient * v1[i];
        v2rho2[i] += term.coefficient * f1[i];
      }
    }
  }

 private:
  struct Term {
    xc_func_type func;
    double coefficient;
    bool gga;
  };
  void release() {
    for (Term& term : terms_) xc_func_end(&term.func);
    terms_.clear();
  }
  std::vector<Term> terms_;
  bool needs_gradient_;
};

}  // namespace

std::unique_ptr<XcEvaluator> make_libxc_evaluator(const Functional& functional) {
  return std::make_unique<LibxcEvaluator>(functional);
}

const char* libxc_version() { return xc_version_string(); }

HybridParameters libxc_hybrid_parameters(const std::string& libxc_name) {
  const int id = xc_functional_get_number(libxc_name.c_str());
  if (id <= 0) throw std::invalid_argument("libxc does not know the functional '" + libxc_name + "'");
  xc_func_type func;
  if (xc_func_init(&func, id, XC_UNPOLARIZED) != 0)
    throw std::runtime_error("libxc could not initialise '" + libxc_name + "'");
  HybridParameters parameters;
  const int family = func.info->family;
  parameters.needs_gradient = family == XC_FAMILY_GGA || family == XC_FAMILY_HYB_GGA;
  if (family == XC_FAMILY_HYB_GGA || family == XC_FAMILY_HYB_LDA) {
    double omega = 0.0, alpha = 0.0, beta = 0.0;
    xc_hyb_cam_coef(&func, &omega, &alpha, &beta);
    parameters.exact_exchange = alpha + beta;
    parameters.long_range_exchange = (beta != 0.0) ? -beta : 0.0;
    parameters.omega = (beta != 0.0) ? omega : 0.0;
  }
  xc_func_end(&func);
  return parameters;
}

}  // namespace ellidft
