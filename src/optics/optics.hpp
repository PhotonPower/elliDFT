#pragma once

// From oligomer polarisabilities to optical constants (specification,
// section 8, "Vom Oligomer zu n und k"). All formulas and conversion factors
// are those given there.

#include <complex>
#include <vector>

#include <Eigen/Dense>

#include "core/molecule.hpp"

namespace ellidft::optics {

// Conversion factors of the specification, section 8.
inline constexpr double angstrom3_per_au = 0.148184711;   // alpha' [A^3] = factor * alpha [a.u.]
inline constexpr double avogadro_factor = 0.602214076;    // N [1/A^3] = factor * rho [g/cm^3] / M [g/mol]
inline constexpr double ev_nm = 1239.841984;              // photon energy [eV] = factor / wavelength [nm]

// Number density of repeat units in 1/Angstrom^3 from the mass density
// (g/cm^3) and the molar mass of the repeat unit (g/mol).
double number_density(double density_g_cm3, double molar_mass_g_mol);

// Photon energy in Hartree for a vacuum wavelength in nm, and back.
double photon_energy_hartree(double wavelength_nm);
double wavelength_nm(double photon_energy_hartree);

// Lorentz-Lorenz relation for a (complex) polarisability volume alpha' in
// Angstrom^3 and a number density N in 1/Angstrom^3:
//
//   (n^2 - 1) / (n^2 + 2) = (4 pi / 3) N alpha'.
//
// Returns the complex refractive index n + i k with n > 0 and k >= 0 for an
// absorbing medium (Im alpha' >= 0). Throws std::domain_error if
// (4 pi / 3) N alpha' is 1 (polarisation catastrophe).
std::complex<double> lorentz_lorenz(std::complex<double> alpha_volume, double number_density);

// Inverse relation: polarisability volume for a given refractive index.
std::complex<double> lorentz_lorenz_inverse(std::complex<double> refractive_index, double number_density);

// Rotation into the chain coordinate system of an oligomer: the columns are
// the principal axes of the distribution of the atoms (all atoms with equal
// weight), ordered so that the third axis (z) is the direction of largest
// extent, i.e. the chain axis. Right-handed.
Eigen::Matrix3d chain_frame(const Molecule& molecule);

// A tensor given in the laboratory frame, expressed in the frame whose axes
// are the columns of `frame`: R^T alpha R.
Eigen::Matrix3cd rotate_into_frame(const Eigen::Matrix3cd& alpha, const Eigen::Matrix3d& frame);

// Polarisability per repeat unit as differences of consecutive oligomers,
// alpha(m) - alpha(m-1). Input: tensors of the oligomers with 1, 2, ... units
// in a common (chain) frame; output: one increment per pair, the first for
// m = 2. The whole series is returned so that the convergence is visible.
std::vector<Eigen::Matrix3cd> increments(const std::vector<Eigen::Matrix3cd>& oligomer_tensors);

// Uniaxial averages of a tensor in the chain frame (z = chain axis).
struct UniaxialPolarizability {
  std::complex<double> parallel;       // alpha_zz
  std::complex<double> perpendicular;  // (alpha_xx + alpha_yy) / 2
  std::complex<double> mean() const { return (parallel + 2.0 * perpendicular) / 3.0; }
  std::complex<double> anisotropy() const { return parallel - perpendicular; }
};
UniaxialPolarizability uniaxial(const Eigen::Matrix3cd& alpha_chain_frame);

// Polarisabilities for the ordinary and extraordinary ray of a uniaxial film
// with order parameter S of the chain axis with respect to the film normal:
//   alpha_e = mean + 2/3 S anisotropy,  alpha_o = mean - 1/3 S anisotropy.
struct OrientedPolarizability {
  std::complex<double> ordinary, extraordinary;
};
OrientedPolarizability oriented(const UniaxialPolarizability& alpha, double order_parameter);

// Optical constants at one wavelength.
struct OpticalConstants {
  double wavelength_nm = 0.0;
  std::complex<double> isotropic;       // n + i k of the isotropic material
  std::complex<double> ordinary;        // uniaxial film; equal to isotropic for S = 0
  std::complex<double> extraordinary;
};

// Optical constants from the polarisability per repeat unit in atomic units
// (chain frame). The local-field correction is applied separately to the
// ordinary and extraordinary component, which is a first approximation for an
// anisotropic medium (specification, section 8); for S = 0 it is exact
// within the Lorentz-Lorenz model.
OpticalConstants optical_constants(double wavelength_nm, const Eigen::Matrix3cd& alpha_unit_au,
                                   double number_density, double order_parameter = 0.0);

}  // namespace ellidft::optics
