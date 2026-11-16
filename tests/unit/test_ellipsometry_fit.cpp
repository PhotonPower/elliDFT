// Module ellipsometry: reading measurements and fitting a layer model.
// The "measurements" are computed from known models, so the fit has a known
// answer.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <cstdint>
#include <memory>
#include <sstream>

#include "ellipsometry/analysis.hpp"

using namespace ellidft::ellipsometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

// Absorbing substrate with some dispersion (invented, silicon-like in size).
std::shared_ptr<Material> substrate() {
  return std::make_shared<TabulatedMaterial>(std::vector<double>{300.0, 500.0, 700.0, 1000.0},
                                             std::vector<Complex>{{5.0, 4.2}, {4.3, 0.07}, {3.8, 0.01}, {3.6, 0.0}});
}

// Film (Cauchy) on a thin oxide (fixed) on the substrate.
struct Sample {
  Stack stack;
  std::shared_ptr<CauchyMaterial> film;
};
Sample make_sample(double thickness, double A, double B) {
  Sample s;
  s.film = std::make_shared<CauchyMaterial>(A, B, 0.0);
  s.stack.ambient = std::make_shared<ConstantMaterial>(1.0);
  s.stack.layers = {{s.film, thickness}, {std::make_shared<CauchyMaterial>(1.45, 0.0036, 0.0), 2.0}};
  s.stack.substrate = substrate();
  return s;
}

std::vector<Measurement> synthetic(const Stack& stack, double noise_deg = 0.0) {
  // deterministic pseudo-random numbers: sum of 12 uniform numbers minus 6 is
  // approximately standard normal
  std::uint64_t state = 20261114;
  auto uniform = [&state] {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return double(state >> 11) / double(1ULL << 53);
  };
  auto normal = [&uniform] {
    double sum = -6.0;
    for (int i = 0; i < 12; ++i) sum += uniform();
    return sum;
  };
  std::vector<Measurement> data;
  for (double angle : {60.0, 70.0})
    for (double wavelength = 400.0; wavelength <= 900.0; wavelength += 20.0) {
      const PsiDelta pd = psi_delta(stack, wavelength, angle);
      Measurement m;
      m.wavelength_nm = wavelength;
      m.angle_deg = angle;
      m.psi_deg = pd.psi_deg + noise_deg * normal();
      m.delta_deg = std::fmod(pd.delta_deg + noise_deg * normal() + 360.0, 360.0);
      data.push_back(m);
    }
  return data;
}

std::vector<FitParameter> film_parameters(Sample& s) {
  std::vector<FitParameter> p = {thickness_parameter(s.stack, 0, "thickness"),
                                 material_parameter(s.film, 0, "A"), material_parameter(s.film, 1, "B")};
  p[1].lower = 1.0;
  p[2].lower = 0.0;
  return p;
}

}  // namespace

TEST_CASE("measurement tables are read", "[ellipsometry]") {
  std::istringstream in("# wavelength angle psi delta\n"
                        "\n"
                        "500 70 12.5 170.25\n"
                        "  600.5, 65; 30 -10   0.02 0.05\n"
                        "# end\n");
  const auto data = read_measurements(in);
  REQUIRE(data.size() == 2);
  CHECK(data[0].wavelength_nm == 500.0);
  CHECK(data[0].angle_deg == 70.0);
  CHECK(data[0].psi_deg == 12.5);
  CHECK(data[0].delta_deg == 170.25);
  CHECK(data[0].sigma_psi_deg == 1.0);
  CHECK(data[1].wavelength_nm == 600.5);
  CHECK(data[1].delta_deg == 350.0);   // brought into [0, 360)
  CHECK(data[1].sigma_psi_deg == 0.02);
  CHECK(data[1].sigma_delta_deg == 0.05);

  std::istringstream shifted("500 70 12.5 170.25\n");
  CHECK(read_measurements(shifted, 180.0).front().delta_deg == 350.25);

  for (const char* bad : {"500 70 12.5\n", "500 70 12.5 abc\n", "500 70 12.5 10 0.1\n", "-500 70 12.5 10\n",
                          "500 95 12.5 10\n", "500 70 120 10\n", "500 70 12.5 10 0 0.1\n", "# nothing\n"}) {
    std::istringstream bad_in(bad);
    CHECK_THROWS_AS(read_measurements(bad_in), std::runtime_error);
  }
  CHECK_THROWS_AS(read_measurements_file("/nonexistent/file.txt"), std::runtime_error);
}

TEST_CASE("fit recovers thickness and dispersion from exact data", "[ellipsometry]") {
  const Sample truth = make_sample(120.0, 1.57, 0.008);
  const auto data = synthetic(truth.stack);

  for (Residuals kind : {Residuals::NCS, Residuals::PsiDelta}) {
    Sample s = make_sample(105.0, 1.50, 0.0);
    const StackFit result = fit_stack(s.stack, film_parameters(s), data, kind);
    CHECK(result.fit.converged);
    CHECK_THAT(result.fit.parameters[0], WithinAbs(120.0, 1e-6));
    CHECK_THAT(result.fit.parameters[1], WithinAbs(1.57, 1e-8));
    CHECK_THAT(result.fit.parameters[2], WithinAbs(0.008, 1e-8));
    // the stack holds the result
    CHECK(s.stack.layers[0].thickness_nm == result.fit.parameters[0]);
    CHECK(s.film->parameter(0) == result.fit.parameters[1]);
    CHECK(result.deviation.rms_psi_deg < 1e-7);
    CHECK(result.deviation.rms_delta_deg < 1e-6);
    CHECK(result.names == std::vector<std::string>{"thickness", "A", "B"});
  }
}

TEST_CASE("fit of noisy data: result within the reported errors", "[ellipsometry]") {
  const Sample truth = make_sample(120.0, 1.57, 0.008);
  const double noise = 0.05;   // degrees
  const auto data = synthetic(truth.stack, noise);

  Sample s = make_sample(110.0, 1.52, 0.004);
  const StackFit result = fit_stack(s.stack, film_parameters(s), data, Residuals::PsiDelta);
  CHECK(result.fit.converged);
  REQUIRE(result.fit.covariance_valid);
  const double truth_values[3] = {120.0, 1.57, 0.008};
  for (int j = 0; j < 3; ++j) {
    INFO(result.names[static_cast<std::size_t>(j)] << " = " << result.fit.parameters[j] << " +- " << result.fit.standard_errors[j]);
    CHECK(result.fit.standard_errors[j] > 0.0);
    CHECK(std::abs(result.fit.parameters[j] - truth_values[j]) < 4.0 * result.fit.standard_errors[j]);
  }
  // unit weights of 1 degree: the residual per degree of freedom estimates the noise
  CHECK_THAT(result.rms_residual, WithinRel(noise, 0.25));
  // the correlation matrix is symmetric with unit diagonal and entries in [-1, 1]
  for (int i = 0; i < 3; ++i) {
    CHECK_THAT(result.fit.correlation(i, i), WithinRel(1.0, 1e-12));
    for (int j = 0; j < i; ++j) {
      CHECK(std::abs(result.fit.correlation(i, j)) < 1.0);
      CHECK_THAT(result.fit.correlation(i, j), WithinAbs(result.fit.correlation(j, i), 1e-12));
    }
  }
}

TEST_CASE("fit with several start values finds a thick film", "[ellipsometry]") {
  // For a thick film the fit has many local minima in the thickness.
  const Sample truth = make_sample(1500.0, 1.57, 0.008);
  const auto data = synthetic(truth.stack);

  Sample s = make_sample(100.0, 1.57, 0.008);
  std::vector<FitParameter> p = {thickness_parameter(s.stack, 0, "thickness")};
  for (double start = 200.0; start <= 2500.0; start += 100.0) p[0].alternative_starts.push_back(start);
  const StackFit result = fit_stack(s.stack, p, data);
  CHECK(result.n_starts == 25);
  CHECK_THAT(result.fit.parameters[0], WithinAbs(1500.0, 1e-5));
  CHECK(result.deviation.rms_ncs < 1e-8);
}

TEST_CASE("fit of the density of a material given by its polarisability", "[ellipsometry]") {
  // Film whose dispersion is fixed by alpha'(lambda); only the density and
  // the thickness are adjusted.
  const std::vector<double> wavelengths = {350.0, 500.0, 700.0, 950.0};
  const std::vector<Complex> alpha = {{14.4, 0.0}, {13.6, 0.0}, {13.2, 0.0}, {13.0, 0.0}};
  auto make = [&](double thickness, double density) {
    Stack stack;
    stack.ambient = std::make_shared<ConstantMaterial>(1.0);
    auto film = std::make_shared<LorentzLorenzMaterial>(wavelengths, alpha, 104.15, density);
    stack.layers = {{film, thickness}};
    stack.substrate = substrate();
    return std::pair{stack, film};
  };
  auto [truth, truth_film] = make(85.0, 1.047);
  const auto data = synthetic(truth);

  auto [stack, film] = make(70.0, 1.2);
  std::vector<FitParameter> p = {thickness_parameter(stack, 0, "thickness"), material_parameter(film, 0, "density")};
  p[1].lower = 0.5;
  p[1].upper = 2.0;
  const StackFit result = fit_stack(stack, p, data);
  CHECK(result.fit.converged);
  CHECK_THAT(result.fit.parameters[0], WithinAbs(85.0, 1e-6));
  CHECK_THAT(result.fit.parameters[1], WithinAbs(1.047, 1e-8));
}

TEST_CASE("Delta near 0 and 360 degrees is compared without a jump", "[ellipsometry]") {
  // bare transparent substrate above the Brewster angle: Delta = 0 (or 360)
  Stack stack{std::make_shared<ConstantMaterial>(1.0), {}, std::make_shared<ConstantMaterial>(1.5)};
  Measurement m;
  m.wavelength_nm = 500.0;
  m.angle_deg = 70.0;
  m.psi_deg = psi_delta(stack, 500.0, 70.0).psi_deg;
  m.delta_deg = 359.9;
  const Eigen::VectorXd r = residual_vector(stack, {m}, Residuals::PsiDelta);
  CHECK_THAT(std::abs(r[1]), WithinAbs(0.1, 1e-9));
  const Deviations d = deviations(stack, {m});
  CHECK_THAT(d.max_delta_deg, WithinAbs(0.1, 1e-9));
  CHECK(residual_vector(stack, {m}, Residuals::NCS).size() == 3);
}

TEST_CASE("fit rejects invalid requests", "[ellipsometry]") {
  Sample s = make_sample(100.0, 1.5, 0.0);
  const auto data = synthetic(s.stack);
  CHECK_THROWS_AS(fit_stack(s.stack, {}, data), std::invalid_argument);
  CHECK_THROWS_AS(fit_stack(s.stack, film_parameters(s), {}), std::invalid_argument);
  CHECK_THROWS_AS(thickness_parameter(s.stack, 5, "x"), std::out_of_range);
  CHECK_THROWS_AS(material_parameter(s.film, 3, "x"), std::out_of_range);
  CHECK_THROWS_AS(material_parameter(substrate(), 0, "x"), std::out_of_range);
}
