#pragma once

// Dispersion models fitted to a computed (or tabulated) refractive index, so
// that the result of the calculation can enter the evaluation of
// measurements as a Cauchy or Sellmeier material (materials.hpp).

#include <vector>

namespace ellidft::ellipsometry {

struct CauchyFit {
  double A = 0.0, B = 0.0, C = 0.0;   // n = A + B / lambda^2 + C / lambda^4, lambda in micrometres
  double rms = 0.0, largest = 0.0;    // deviation of the model from the data
};
// Linear least squares. With fewer than four points C is fixed at zero.
// Throws std::invalid_argument for fewer than two points.
CauchyFit fit_cauchy(const std::vector<double>& wavelengths_nm, const std::vector<double>& n);

struct SellmeierFit {
  double B = 0.0, C_um2 = 0.0;        // n^2 = 1 + B lambda^2 / (lambda^2 - C), lambda in micrometres
  double rms = 0.0, largest = 0.0;
  bool converged = false;
};
// One term. Start values from the straight line 1 / (n^2 - 1) = 1/B - (C/B) / lambda^2,
// then least squares in n. Needs at least three points with n > 1.
SellmeierFit fit_sellmeier(const std::vector<double>& wavelengths_nm, const std::vector<double>& n);

}  // namespace ellidft::ellipsometry
