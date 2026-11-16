#include "ellipsometry/materials.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "optics/optics.hpp"

namespace ellidft::ellipsometry {

namespace {

// Linear interpolation in a table with ascending abscissae; no extrapolation.
Complex interpolate(const std::vector<double>& x, const std::vector<Complex>& y, double at) {
  if (at < x.front() || at > x.back())
    throw std::out_of_range("material table covers " + std::to_string(x.front()) + " to " +
                            std::to_string(x.back()) + " nm, requested " + std::to_string(at) + " nm");
  const auto upper = std::upper_bound(x.begin(), x.end(), at);
  if (upper == x.end()) return y.back();
  const auto i = static_cast<std::size_t>(upper - x.begin());
  if (i == 0) return y.front();
  const double w = (at - x[i - 1]) / (x[i] - x[i - 1]);
  return (1.0 - w) * y[i - 1] + w * y[i];
}

void check_table(const std::vector<double>& x, std::size_t n_values) {
  if (x.empty() || x.size() != n_values) throw std::invalid_argument("material table: sizes differ or table is empty");
  for (std::size_t i = 1; i < x.size(); ++i)
    if (!(x[i] > x[i - 1])) throw std::invalid_argument("material table: wavelengths must be strictly ascending");
}

}  // namespace

double Material::parameter(std::size_t) const { throw std::out_of_range("material has no parameters"); }
void Material::set_parameter(std::size_t, double) { throw std::out_of_range("material has no parameters"); }

double ConstantMaterial::parameter(std::size_t i) const {
  if (i > 1) throw std::out_of_range("ConstantMaterial: parameter index");
  return i == 0 ? n_ : k_;
}
void ConstantMaterial::set_parameter(std::size_t i, double value) {
  if (i > 1) throw std::out_of_range("ConstantMaterial: parameter index");
  (i == 0 ? n_ : k_) = value;
}

Complex CauchyMaterial::index(double wavelength_nm) const {
  const double l2 = (wavelength_nm / 1000.0) * (wavelength_nm / 1000.0);
  return {a_[0] + a_[1] / l2 + a_[2] / (l2 * l2), 0.0};
}

SellmeierMaterial::SellmeierMaterial(std::vector<double> B, std::vector<double> C_um2)
    : B_(std::move(B)), C_(std::move(C_um2)) {
  if (B_.empty() || B_.size() != C_.size())
    throw std::invalid_argument("SellmeierMaterial: need equally many B and C coefficients");
}

Complex SellmeierMaterial::index(double wavelength_nm) const {
  const double l2 = (wavelength_nm / 1000.0) * (wavelength_nm / 1000.0);
  double n2 = 1.0;
  for (std::size_t i = 0; i < B_.size(); ++i) n2 += B_[i] * l2 / (l2 - C_[i]);
  if (!(n2 > 0.0)) throw std::domain_error("SellmeierMaterial: n^2 is not positive at this wavelength");
  return {std::sqrt(n2), 0.0};
}

std::vector<std::string> SellmeierMaterial::parameter_names() const {
  std::vector<std::string> names;
  for (std::size_t i = 0; i < B_.size(); ++i) {
    names.push_back("B" + std::to_string(i + 1));
    names.push_back("C" + std::to_string(i + 1));
  }
  return names;
}
double SellmeierMaterial::parameter(std::size_t i) const { return (i % 2 == 0 ? B_ : C_).at(i / 2); }
void SellmeierMaterial::set_parameter(std::size_t i, double value) { (i % 2 == 0 ? B_ : C_).at(i / 2) = value; }

TabulatedMaterial::TabulatedMaterial(std::vector<double> wavelengths_nm, std::vector<Complex> indices)
    : wavelengths_(std::move(wavelengths_nm)), indices_(std::move(indices)) {
  check_table(wavelengths_, indices_.size());
}

Complex TabulatedMaterial::index(double wavelength_nm) const {
  return interpolate(wavelengths_, indices_, wavelength_nm);
}

LorentzLorenzMaterial::LorentzLorenzMaterial(std::vector<double> wavelengths_nm, std::vector<Complex> alpha_volume_A3,
                                             double molar_mass_g_mol, double density_g_cm3)
    : wavelengths_(std::move(wavelengths_nm)), alpha_(std::move(alpha_volume_A3)), molar_mass_(molar_mass_g_mol),
      density_(density_g_cm3) {
  check_table(wavelengths_, alpha_.size());
  if (!(molar_mass_ > 0.0)) throw std::invalid_argument("LorentzLorenzMaterial: molar mass must be positive");
}

Complex LorentzLorenzMaterial::index(double wavelength_nm) const {
  const Complex alpha = interpolate(wavelengths_, alpha_, wavelength_nm);
  return optics::lorentz_lorenz(alpha, optics::number_density(density_, molar_mass_));
}

double LorentzLorenzMaterial::parameter(std::size_t i) const {
  if (i != 0) throw std::out_of_range("LorentzLorenzMaterial: parameter index");
  return density_;
}
void LorentzLorenzMaterial::set_parameter(std::size_t i, double value) {
  if (i != 0) throw std::out_of_range("LorentzLorenzMaterial: parameter index");
  density_ = value;
}

namespace {

// Square root of a permittivity with non-negative imaginary part of the index.
Complex index_of(Complex eps) {
  Complex N = std::sqrt(eps);
  if (N.imag() < 0.0 || (N.imag() == 0.0 && N.real() < 0.0)) N = -N;
  return N;
}

}  // namespace

double photon_energy_ev(double wavelength_nm) {
  if (!(wavelength_nm > 0.0)) throw std::invalid_argument("photon_energy_ev: wavelength must be positive");
  return optics::ev_nm / wavelength_nm;
}

LorentzMaterial::LorentzMaterial(double eps_inf, std::vector<Oscillator> oscillators)
    : eps_inf_(eps_inf), oscillators_(std::move(oscillators)) {}

Complex LorentzMaterial::permittivity(double energy_ev) const {
  Complex eps(eps_inf_, 0.0);
  for (const Oscillator& o : oscillators_)
    eps += o.strength_ev2 / Complex(o.energy_ev * o.energy_ev - energy_ev * energy_ev, -o.width_ev * energy_ev);
  return eps;
}

Complex LorentzMaterial::index(double wavelength_nm) const { return index_of(permittivity(photon_energy_ev(wavelength_nm))); }

std::vector<std::string> LorentzMaterial::parameter_names() const {
  std::vector<std::string> names = {"eps_inf"};
  for (std::size_t j = 0; j < oscillators_.size(); ++j) {
    const std::string k = std::to_string(j + 1);
    names.push_back("strength" + k);
    names.push_back("energy" + k);
    names.push_back("width" + k);
  }
  return names;
}

double LorentzMaterial::parameter(std::size_t i) const {
  if (i == 0) return eps_inf_;
  const Oscillator& o = oscillators_.at((i - 1) / 3);
  return (i - 1) % 3 == 0 ? o.strength_ev2 : ((i - 1) % 3 == 1 ? o.energy_ev : o.width_ev);
}

void LorentzMaterial::set_parameter(std::size_t i, double value) {
  if (i == 0) {
    eps_inf_ = value;
    return;
  }
  Oscillator& o = oscillators_.at((i - 1) / 3);
  ((i - 1) % 3 == 0 ? o.strength_ev2 : ((i - 1) % 3 == 1 ? o.energy_ev : o.width_ev)) = value;
}

double kramers_kronig_real_part(const std::function<double(double)>& eps2, double lower, double scale, double energy) {
  if (!(scale > 0.0) || lower < 0.0 || energy < 0.0) throw std::invalid_argument("kramers_kronig_real_part: invalid arguments");
  // The singularity at xi = E is taken out:
  //   P int f(xi) / (xi^2 - E^2) = int [f(xi) - f(E)] / (xi^2 - E^2) + f(E) P int d xi / (xi^2 - E^2),
  //   f(xi) = xi eps2(xi),  P int_a^inf d xi / (xi^2 - E^2) = ln|(a + E) / (a - E)| / (2 E).
  // The first integrand is smooth. xi = lower + scale u / (1 - u), u in [0, 1), composite Simpson rule.
  const bool singular = energy > lower;
  const double fE = singular ? energy * eps2(energy) : 0.0;
  auto integrand = [&](double u) {
    const double xi = lower + scale * u / (1.0 - u);
    const double jacobian = scale / ((1.0 - u) * (1.0 - u));
    const double denominator = xi * xi - energy * energy;
    if (std::abs(denominator) < 1e-10 * (energy * energy + scale * scale)) {
      const double h = 1e-4 * scale;
      // E = 0 and xi -> 0: the integrand is eps2(xi) / xi
      if (!(energy > 1e-6 * scale)) return jacobian * eps2(h) / h;
      // limit of the difference quotient: d(xi eps2) / d xi / (2 E)
      return jacobian * ((energy + h) * eps2(energy + h) - (energy - h) * eps2(energy - h)) / (2.0 * h) / (2.0 * energy);
    }
    return jacobian * (xi * eps2(xi) - fE) / denominator;
  };
  // At u = 1 (xi = infinity) the integrand tends to (lim xi eps2 - f(E)) / scale = -f(E) / scale,
  // because eps2 falls off faster than 1 / xi.
  const int n = 4000;
  double sum = integrand(0.0) - fE / scale;
  for (int k = 1; k < n; ++k) sum += (k % 2 == 1 ? 4.0 : 2.0) * integrand(static_cast<double>(k) / n);
  double value = sum / (3.0 * n);
  if (singular && lower != energy) value += fE * std::log(std::abs((lower + energy) / (lower - energy))) / (2.0 * energy);
  return 2.0 / M_PI * value;
}

TaucLorentzMaterial::TaucLorentzMaterial(double eps_inf, double A, double E0, double C, double Eg) : p_{eps_inf, A, E0, C, Eg} {}

double TaucLorentzMaterial::eps2(double E) const {
  const double A = p_[1], E0 = p_[2], C = p_[3], Eg = p_[4];
  if (!(E > Eg) || !(E > 0.0)) return 0.0;
  const double d = E * E - E0 * E0;
  return A * E0 * C * (E - Eg) * (E - Eg) / (E * (d * d + C * C * E * E));
}

Complex TaucLorentzMaterial::permittivity(double energy_ev) const {
  if (!(p_[2] > 0.0) || !(p_[3] > 0.0) || p_[4] < 0.0)
    throw std::domain_error("TaucLorentzMaterial: E0 and C must be positive, Eg not negative");
  const double eps1 = p_[0] + kramers_kronig_real_part([this](double xi) { return eps2(xi); }, p_[4], p_[2], energy_ev);
  return {eps1, eps2(energy_ev)};
}

Complex TaucLorentzMaterial::index(double wavelength_nm) const { return index_of(permittivity(photon_energy_ev(wavelength_nm))); }

BruggemanMaterial::BruggemanMaterial(std::shared_ptr<Material> a, std::shared_ptr<Material> b, double fraction_b)
    : a_(std::move(a)), b_(std::move(b)), fraction_b_(fraction_b) {
  if (!a_ || !b_) throw std::invalid_argument("BruggemanMaterial: two components are required");
}

Complex BruggemanMaterial::mix(Complex eps_a, Complex eps_b, double fraction_b) {
  if (fraction_b < 0.0 || fraction_b > 1.0) throw std::domain_error("BruggemanMaterial: the fraction must be in [0, 1]");
  const double fa = 1.0 - fraction_b, fb = fraction_b;
  const Complex H = (3.0 * fa - 1.0) * eps_a + (3.0 * fb - 1.0) * eps_b;
  const Complex root = std::sqrt(H * H + 8.0 * eps_a * eps_b);
  const Complex plus = (H + root) / 4.0, minus = (H - root) / 4.0;
  // physical root: absorbing or transparent, never amplifying; for transparent components the positive one
  auto acceptable = [](Complex e) { return e.imag() >= -1e-14 * std::abs(e) && (e.imag() > 1e-14 * std::abs(e) || e.real() > 0.0); };
  if (acceptable(plus) && !acceptable(minus)) return plus;
  if (acceptable(minus) && !acceptable(plus)) return minus;
  return plus.imag() >= minus.imag() ? plus : minus;
}

Complex BruggemanMaterial::index(double wavelength_nm) const {
  const Complex Na = a_->index(wavelength_nm), Nb = b_->index(wavelength_nm);
  return index_of(mix(Na * Na, Nb * Nb, fraction_b_));
}

double BruggemanMaterial::parameter(std::size_t i) const {
  if (i != 0) throw std::out_of_range("BruggemanMaterial: parameter index");
  return fraction_b_;
}
void BruggemanMaterial::set_parameter(std::size_t i, double value) {
  if (i != 0) throw std::out_of_range("BruggemanMaterial: parameter index");
  fraction_b_ = value;
}

}  // namespace ellidft::ellipsometry
