#pragma once

// Evaluation of ellipsometry measurements: measured Psi and Delta, a layer
// model with adjustable parameters, and the fit of the one to the other.

#include <functional>
#include <istream>
#include <limits>
#include <string>
#include <vector>

#include "ellipsometry/least_squares.hpp"
#include "ellipsometry/stack.hpp"

namespace ellidft::ellipsometry {

struct Measurement {
  double wavelength_nm = 0.0;
  double angle_deg = 0.0;
  double psi_deg = 0.0;
  double delta_deg = 0.0;   // in [0, 360), convention of psi_delta() in stack.hpp
  double sigma_psi_deg = 1.0, sigma_delta_deg = 1.0;   // only used with Residuals::PsiDelta
};

// Text table, one measurement per line:
//   wavelength_nm  angle_deg  psi_deg  delta_deg  [sigma_psi_deg  sigma_delta_deg]
// Blank lines and lines starting with '#' are skipped. `delta_offset_deg` is
// added to every Delta (use 180 for data in the convention with the opposite
// sign of r_p); Delta is then brought into [0, 360).
std::vector<Measurement> read_measurements(std::istream& in, double delta_offset_deg = 0.0);
std::vector<Measurement> read_measurements_file(const std::string& path, double delta_offset_deg = 0.0);

// What is compared.
//   NCS:      N = cos 2Psi, C = sin 2Psi cos Delta, S = sin 2Psi sin Delta,
//             three residuals per point, unweighted. Continuous everywhere,
//             also where Delta jumps by 360 degrees or is undefined (Psi -> 0).
//   PsiDelta: (Psi_meas - Psi_calc) / sigma_Psi and the Delta difference,
//             taken into (-180, 180], divided by sigma_Delta; two residuals
//             per point.
enum class Residuals { NCS, PsiDelta };

// An adjustable quantity of the model.
struct FitParameter {
  std::string name;
  std::function<double()> get;
  std::function<void(double)> set;
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();
  // Start values to try in addition to the current one (the fit is run from
  // every combination and the best result kept). Meant for film thicknesses,
  // where the fit has several minima.
  std::vector<double> alternative_starts;
};
FitParameter thickness_parameter(Stack& stack, std::size_t layer, std::string name);
FitParameter material_parameter(const std::shared_ptr<Material>& material, std::size_t index, std::string name);

std::vector<PsiDelta> simulate(const Stack& stack, const std::vector<Measurement>& points);
Eigen::VectorXd residual_vector(const Stack& stack, const std::vector<Measurement>& data, Residuals kind);

struct Deviations {
  double rms_psi_deg = 0.0, rms_delta_deg = 0.0, rms_ncs = 0.0;
  double max_psi_deg = 0.0, max_delta_deg = 0.0;
};
Deviations deviations(const Stack& stack, const std::vector<Measurement>& data);

struct StackFit {
  LeastSquaresResult fit;              // best of all starts
  std::vector<std::string> names;
  int n_starts = 0;
  // sqrt(S / (n - m)): root mean square of the residuals per degree of freedom.
  double rms_residual = 0.0;
  Deviations deviation;
};

// Adjusts the parameters; on return the stack holds the fitted values.
StackFit fit_stack(Stack& stack, std::vector<FitParameter> parameters, const std::vector<Measurement>& data,
                   Residuals kind = Residuals::NCS, const LeastSquaresOptions& options = {});

}  // namespace ellidft::ellipsometry
