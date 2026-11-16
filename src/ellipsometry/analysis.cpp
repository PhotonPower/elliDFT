#include "ellipsometry/analysis.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ellidft::ellipsometry {

namespace {

constexpr double degree = M_PI / 180.0;
constexpr std::size_t max_starts = 10000;

double wrap_360(double angle) {
  angle = std::fmod(angle, 360.0);
  return angle < 0.0 ? angle + 360.0 : angle;
}

// Difference of two angles in (-180, 180].
double angle_difference(double a, double b) {
  double d = std::fmod(a - b, 360.0);
  if (d > 180.0) d -= 360.0;
  if (d <= -180.0) d += 360.0;
  return d;
}

struct NCS {
  double n, c, s;
};
NCS ncs(double psi_deg, double delta_deg) {
  const double two_psi = 2.0 * psi_deg * degree, delta = delta_deg * degree;
  return {std::cos(two_psi), std::sin(two_psi) * std::cos(delta), std::sin(two_psi) * std::sin(delta)};
}

}  // namespace

std::vector<Measurement> read_measurements(std::istream& in, double delta_offset_deg) {
  std::vector<Measurement> data;
  std::string line;
  int line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    for (char& c : line)
      if (c == ',' || c == ';') c = ' ';
    std::istringstream fields(line);
    std::vector<double> values;
    for (std::string token; fields >> token;) {
      std::size_t used = 0;
      double v = 0.0;
      try {
        v = std::stod(token, &used);
      } catch (const std::exception&) {
        used = 0;
      }
      if (used != token.size()) throw std::runtime_error("measurement table, line " + std::to_string(line_number) + ": '" + token + "' is not a number");
      values.push_back(v);
    }
    if (values.size() != 4 && values.size() != 6)
      throw std::runtime_error("measurement table, line " + std::to_string(line_number) + ": expected 4 or 6 columns, found " + std::to_string(values.size()));
    Measurement m;
    m.wavelength_nm = values[0];
    m.angle_deg = values[1];
    m.psi_deg = values[2];
    m.delta_deg = wrap_360(values[3] + delta_offset_deg);
    if (values.size() == 6) {
      m.sigma_psi_deg = values[4];
      m.sigma_delta_deg = values[5];
    }
    const std::string where = "measurement table, line " + std::to_string(line_number) + ": ";
    if (!(m.wavelength_nm > 0.0)) throw std::runtime_error(where + "wavelength must be positive");
    if (!(m.angle_deg >= 0.0 && m.angle_deg < 90.0)) throw std::runtime_error(where + "angle of incidence must be in [0, 90)");
    if (!(m.psi_deg >= 0.0 && m.psi_deg <= 90.0)) throw std::runtime_error(where + "Psi must be in [0, 90]");
    if (!(m.sigma_psi_deg > 0.0 && m.sigma_delta_deg > 0.0)) throw std::runtime_error(where + "standard deviations must be positive");
    data.push_back(m);
  }
  if (data.empty()) throw std::runtime_error("measurement table: no data");
  return data;
}

std::vector<Measurement> read_measurements_file(const std::string& path, double delta_offset_deg) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  try {
    return read_measurements(in, delta_offset_deg);
  } catch (const std::exception& e) {
    throw std::runtime_error(path + ": " + e.what());
  }
}

FitParameter thickness_parameter(Stack& stack, std::size_t layer, std::string name) {
  if (layer >= stack.layers.size()) throw std::out_of_range("thickness_parameter: no such layer");
  FitParameter p;
  p.name = std::move(name);
  p.get = [&stack, layer] { return stack.layers[layer].thickness_nm; };
  p.set = [&stack, layer](double value) { stack.layers[layer].thickness_nm = value; };
  p.lower = 0.0;
  return p;
}

FitParameter material_parameter(const std::shared_ptr<Material>& material, std::size_t index, std::string name) {
  if (!material || index >= material->parameter_names().size())
    throw std::out_of_range("material_parameter: no such parameter");
  FitParameter p;
  p.name = std::move(name);
  p.get = [material, index] { return material->parameter(index); };
  p.set = [material, index](double value) { material->set_parameter(index, value); };
  return p;
}

std::vector<PsiDelta> simulate(const Stack& stack, const std::vector<Measurement>& points) {
  std::vector<PsiDelta> result;
  result.reserve(points.size());
  for (const Measurement& m : points) result.push_back(psi_delta(stack, m.wavelength_nm, m.angle_deg));
  return result;
}

Eigen::VectorXd residual_vector(const Stack& stack, const std::vector<Measurement>& data, Residuals kind) {
  const Eigen::Index per_point = kind == Residuals::NCS ? 3 : 2;
  Eigen::VectorXd r(per_point * static_cast<Eigen::Index>(data.size()));
  Eigen::Index at = 0;
  for (const Measurement& m : data) {
    const PsiDelta c = psi_delta(stack, m.wavelength_nm, m.angle_deg);
    if (kind == Residuals::NCS) {
      const NCS measured = ncs(m.psi_deg, m.delta_deg), calculated = ncs(c.psi_deg, c.delta_deg);
      r[at++] = measured.n - calculated.n;
      r[at++] = measured.c - calculated.c;
      r[at++] = measured.s - calculated.s;
    } else {
      r[at++] = (m.psi_deg - c.psi_deg) / m.sigma_psi_deg;
      r[at++] = angle_difference(m.delta_deg, c.delta_deg) / m.sigma_delta_deg;
    }
  }
  return r;
}

Deviations deviations(const Stack& stack, const std::vector<Measurement>& data) {
  Deviations d;
  double sum_psi = 0.0, sum_delta = 0.0, sum_ncs = 0.0;
  for (const Measurement& m : data) {
    const PsiDelta c = psi_delta(stack, m.wavelength_nm, m.angle_deg);
    const double dpsi = m.psi_deg - c.psi_deg, ddelta = angle_difference(m.delta_deg, c.delta_deg);
    const NCS a = ncs(m.psi_deg, m.delta_deg), b = ncs(c.psi_deg, c.delta_deg);
    sum_psi += dpsi * dpsi;
    sum_delta += ddelta * ddelta;
    sum_ncs += (a.n - b.n) * (a.n - b.n) + (a.c - b.c) * (a.c - b.c) + (a.s - b.s) * (a.s - b.s);
    d.max_psi_deg = std::max(d.max_psi_deg, std::abs(dpsi));
    d.max_delta_deg = std::max(d.max_delta_deg, std::abs(ddelta));
  }
  const double n = static_cast<double>(data.size());
  if (n > 0.0) {
    d.rms_psi_deg = std::sqrt(sum_psi / n);
    d.rms_delta_deg = std::sqrt(sum_delta / n);
    d.rms_ncs = std::sqrt(sum_ncs / (3.0 * n));
  }
  return d;
}

StackFit fit_stack(Stack& stack, std::vector<FitParameter> parameters, const std::vector<Measurement>& data,
                   Residuals kind, const LeastSquaresOptions& options) {
  if (parameters.empty()) throw std::invalid_argument("fit_stack: no parameters to adjust");
  if (data.empty()) throw std::invalid_argument("fit_stack: no measurements");
  const Eigen::Index m = static_cast<Eigen::Index>(parameters.size());
  Eigen::VectorXd lower(m), upper(m);
  std::vector<std::vector<double>> starts(parameters.size());
  std::size_t n_starts = 1;
  for (Eigen::Index j = 0; j < m; ++j) {
    const FitParameter& p = parameters[static_cast<std::size_t>(j)];
    lower[j] = p.lower;
    upper[j] = p.upper;
    auto& s = starts[static_cast<std::size_t>(j)];
    s.push_back(p.get());
    s.insert(s.end(), p.alternative_starts.begin(), p.alternative_starts.end());
    n_starts *= s.size();
    if (n_starts > max_starts) throw std::invalid_argument("fit_stack: more than 10000 combinations of start values");
  }

  auto apply = [&](const Eigen::VectorXd& values) {
    for (Eigen::Index j = 0; j < m; ++j) parameters[static_cast<std::size_t>(j)].set(values[j]);
  };
  const ResidualFunction residuals = [&](const Eigen::VectorXd& values) {
    apply(values);
    return residual_vector(stack, data, kind);
  };

  StackFit result;
  for (const FitParameter& p : parameters) result.names.push_back(p.name);
  result.n_starts = static_cast<int>(n_starts);
  bool have_result = false;
  std::string last_error;
  for (std::size_t combination = 0; combination < n_starts; ++combination) {
    Eigen::VectorXd start(m);
    std::size_t rest = combination;
    for (Eigen::Index j = 0; j < m; ++j) {
      const auto& s = starts[static_cast<std::size_t>(j)];
      start[j] = s[rest % s.size()];
      rest /= s.size();
    }
    try {
      LeastSquaresResult fit = levenberg_marquardt(residuals, start, lower, upper, options);
      if (!have_result || fit.sum_of_squares < result.fit.sum_of_squares) {
        result.fit = std::move(fit);
        have_result = true;
      }
    } catch (const std::exception& e) {
      // A start outside the validity of a material model (for example a
      // Sellmeier pole inside the spectrum) is skipped; the others count.
      last_error = e.what();
    }
  }
  if (!have_result) throw std::runtime_error("fit_stack: no start value could be evaluated: " + last_error);

  apply(result.fit.parameters);
  const Eigen::Index n = residual_vector(stack, data, kind).size();
  result.rms_residual = n > m ? std::sqrt(result.fit.sum_of_squares / static_cast<double>(n - m)) : 0.0;
  result.deviation = deviations(stack, data);
  return result;
}

}  // namespace ellidft::ellipsometry
