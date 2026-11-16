#pragma once

// Optical materials for the ellipsometry model: complex refractive index
// N = n + i k as a function of the vacuum wavelength. Convention: time
// dependence exp(-i omega t), so that k >= 0 means absorption.

#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace ellidft::ellipsometry {

using Complex = std::complex<double>;

class Material {
 public:
  virtual ~Material() = default;

  // Complex refractive index at a vacuum wavelength in nm.
  virtual Complex index(double wavelength_nm) const = 0;

  // Adjustable parameters of the model (for fits): names, values, setter.
  virtual std::vector<std::string> parameter_names() const { return {}; }
  virtual double parameter(std::size_t) const;
  virtual void set_parameter(std::size_t, double);
};

// Wavelength-independent index.  Parameters: n, k.
class ConstantMaterial final : public Material {
 public:
  ConstantMaterial(double n, double k = 0.0) : n_(n), k_(k) {}
  Complex index(double) const override { return {n_, k_}; }
  std::vector<std::string> parameter_names() const override { return {"n", "k"}; }
  double parameter(std::size_t i) const override;
  void set_parameter(std::size_t i, double value) override;

 private:
  double n_, k_;
};

// Cauchy model of a transparent material, the usual empirical form
//   n(lambda) = A + B / lambda^2 + C / lambda^4,   lambda in micrometres, k = 0.
// Parameters: A, B, C.
class CauchyMaterial final : public Material {
 public:
  CauchyMaterial(double A, double B, double C = 0.0) : a_{A, B, C} {}
  Complex index(double wavelength_nm) const override;
  std::vector<std::string> parameter_names() const override { return {"A", "B", "C"}; }
  double parameter(std::size_t i) const override { return a_.at(i); }
  void set_parameter(std::size_t i, double value) override { a_.at(i) = value; }

 private:
  std::vector<double> a_;
};

// Sellmeier model of a transparent material,
//   n^2 = 1 + sum_i B_i lambda^2 / (lambda^2 - C_i),   lambda in micrometres, C_i in um^2
// (form of the fits in the refractiveindex.info database, see
// docs/sources/README.md). Parameters: B1, C1, B2, C2, ...
class SellmeierMaterial final : public Material {
 public:
  SellmeierMaterial(std::vector<double> B, std::vector<double> C_um2);
  Complex index(double wavelength_nm) const override;
  std::vector<std::string> parameter_names() const override;
  double parameter(std::size_t i) const override;
  void set_parameter(std::size_t i, double value) override;

 private:
  std::vector<double> B_, C_;
};

// Tabulated n and k, linear interpolation in the wavelength. Outside the
// table the request is an error (no extrapolation).
class TabulatedMaterial final : public Material {
 public:
  TabulatedMaterial(std::vector<double> wavelengths_nm, std::vector<Complex> indices);
  Complex index(double wavelength_nm) const override;

 private:
  std::vector<double> wavelengths_;
  std::vector<Complex> indices_;
};

// Material defined by the polarisability volume per repeat unit alpha'(lambda)
// (Angstrom^3, complex), as computed by elliDFT, and the Lorentz-Lorenz
// relation at a given mass density. Parameter: density (g/cm^3). This links
// the fit of a measurement to the first-principles polarisability: the
// dispersion is fixed by the calculation, the density is adjusted.
class LorentzLorenzMaterial final : public Material {
 public:
  LorentzLorenzMaterial(std::vector<double> wavelengths_nm, std::vector<Complex> alpha_volume_A3,
                        double molar_mass_g_mol, double density_g_cm3);
  Complex index(double wavelength_nm) const override;
  std::vector<std::string> parameter_names() const override { return {"density"}; }
  double parameter(std::size_t i) const override;
  void set_parameter(std::size_t i, double value) override;

 private:
  std::vector<double> wavelengths_;
  std::vector<Complex> alpha_;
  double molar_mass_, density_;
};

// Photon energy in eV for a vacuum wavelength in nm: 1239.841984 / lambda
// (specification, section 8).
double photon_energy_ev(double wavelength_nm);

// Sum of Lorentz oscillators (Wikipedia, "Lorentz oscillator model", time
// dependence exp(-i omega t)), written in photon energies:
//   eps(E) = eps_inf + sum_j s_j / (E_j^2 - E^2 - i Gamma_j E),
// s_j in eV^2, E_j and Gamma_j in eV. N = sqrt(eps) with k >= 0.
// Parameters: eps_inf, then s1, E1, Gamma1, s2, ...
class LorentzMaterial final : public Material {
 public:
  struct Oscillator {
    double strength_ev2, energy_ev, width_ev;
  };
  LorentzMaterial(double eps_inf, std::vector<Oscillator> oscillators);
  Complex permittivity(double energy_ev) const;
  Complex index(double wavelength_nm) const override;
  std::vector<std::string> parameter_names() const override;
  double parameter(std::size_t i) const override;
  void set_parameter(std::size_t i, double value) override;

 private:
  double eps_inf_;
  std::vector<Oscillator> oscillators_;
};

// Real part of a permittivity from its imaginary part by the Kramers-Kronig
// relation,
//   eps1(E) - eps_inf = (2 / pi) P integral_lower^infinity xi eps2(xi) / (xi^2 - E^2) d xi,
// by numerical quadrature. eps2 must vanish below `lower`, be continuous and
// fall off faster than 1 / xi for large xi; `scale` is a typical energy of
// its structure (position of the absorption).
double kramers_kronig_real_part(const std::function<double(double)>& eps2, double lower, double scale, double energy);

// Tauc-Lorentz model of an amorphous absorber (Wikipedia, "Tauc-Lorentz
// model", there after Jellison and Modine 1996):
//   eps2(E) = A E0 C (E - Eg)^2 / (E [(E^2 - E0^2)^2 + C^2 E^2])   for E > Eg, else 0;
//   eps1 from the Kramers-Kronig relation above with eps_inf.
// Energies in eV. The real part is integrated numerically; the closed form
// of the same source is used as a check in the tests only.
// Parameters: eps_inf, A, E0, C, Eg.
class TaucLorentzMaterial final : public Material {
 public:
  TaucLorentzMaterial(double eps_inf, double A, double E0, double C, double Eg);
  double eps2(double energy_ev) const;
  Complex permittivity(double energy_ev) const;
  Complex index(double wavelength_nm) const override;
  std::vector<std::string> parameter_names() const override { return {"eps_inf", "A", "E0", "C", "Eg"}; }
  double parameter(std::size_t i) const override { return p_.at(i); }
  void set_parameter(std::size_t i, double value) override { p_.at(i) = value; }

 private:
  std::vector<double> p_;
};

// Mixture of two materials in the effective-medium approximation of
// Bruggeman for spherical inclusions (Wikipedia, "Effective medium
// approximations"):
//   eps = [H + sqrt(H^2 + 8 eps_a eps_b)] / 4,  H = (3 f_a - 1) eps_a + (3 f_b - 1) eps_b,
// f_a + f_b = 1. Of the two roots the one with Im eps >= 0 (and positive real
// part for transparent components) is taken; the source notes that the sign
// of the root must be chosen for the correct imaginary part.
// Uses: surface roughness (film material mixed with the ambient), porous or
// intermixed layers. Parameter: fraction_b. The parameters of the two
// components belong to those materials.
class BruggemanMaterial final : public Material {
 public:
  BruggemanMaterial(std::shared_ptr<Material> a, std::shared_ptr<Material> b, double fraction_b);
  static Complex mix(Complex eps_a, Complex eps_b, double fraction_b);
  Complex index(double wavelength_nm) const override;
  std::vector<std::string> parameter_names() const override { return {"fraction_b"}; }
  double parameter(std::size_t i) const override;
  void set_parameter(std::size_t i, double value) override;

 private:
  std::shared_ptr<Material> a_, b_;
  double fraction_b_;
};

}  // namespace ellidft::ellipsometry
