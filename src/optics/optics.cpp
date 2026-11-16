#include "optics/optics.hpp"

#include <cmath>
#include <stdexcept>

#include "core/units.hpp"

namespace ellidft::optics {

double number_density(double density_g_cm3, double molar_mass_g_mol) {
  if (!(density_g_cm3 > 0.0) || !(molar_mass_g_mol > 0.0))
    throw std::invalid_argument("number_density: density and molar mass must be positive");
  return avogadro_factor * density_g_cm3 / molar_mass_g_mol;
}

double photon_energy_hartree(double wavelength) {
  if (!(wavelength > 0.0)) throw std::invalid_argument("photon_energy_hartree: wavelength must be positive");
  return ev_nm / wavelength / units::hartree_in_ev;
}

double wavelength_nm(double photon_energy) {
  if (!(photon_energy > 0.0)) throw std::invalid_argument("wavelength_nm: photon energy must be positive");
  return ev_nm / (photon_energy * units::hartree_in_ev);
}

std::complex<double> lorentz_lorenz(std::complex<double> alpha_volume, double density) {
  // (n^2 - 1)/(n^2 + 2) = x  =>  n^2 = (1 + 2x)/(1 - x),  x = (4 pi / 3) N alpha'.
  const std::complex<double> x = 4.0 * M_PI / 3.0 * density * alpha_volume;
  if (std::abs(1.0 - x) < 1e-12) throw std::domain_error("lorentz_lorenz: (4 pi / 3) N alpha' equals 1");
  const std::complex<double> n2 = (1.0 + 2.0 * x) / (1.0 - x);
  std::complex<double> n = std::sqrt(n2);  // principal root: real part >= 0
  if (n.real() < 0.0 || (n.real() == 0.0 && n.imag() < 0.0)) n = -n;
  return n;
}

std::complex<double> lorentz_lorenz_inverse(std::complex<double> refractive_index, double density) {
  const std::complex<double> n2 = refractive_index * refractive_index;
  return (n2 - 1.0) / (n2 + 2.0) * 3.0 / (4.0 * M_PI * density);
}

Eigen::Matrix3d chain_frame(const Molecule& molecule) {
  if (molecule.n_atoms() < 2) throw std::invalid_argument("chain_frame: at least two atoms are needed");
  Eigen::Vector3d center = Eigen::Vector3d::Zero();
  for (const Atom& atom : molecule.atoms) center += Eigen::Vector3d(atom.r[0], atom.r[1], atom.r[2]);
  center /= static_cast<double>(molecule.n_atoms());
  Eigen::Matrix3d spread = Eigen::Matrix3d::Zero();
  for (const Atom& atom : molecule.atoms) {
    const Eigen::Vector3d d = Eigen::Vector3d(atom.r[0], atom.r[1], atom.r[2]) - center;
    spread += d * d.transpose();
  }
  // Eigenvalues ascending: the last eigenvector is the direction of largest extent.
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(spread);
  Eigen::Matrix3d frame = eigen.eigenvectors();
  if (frame.determinant() < 0.0) frame.col(0) = -frame.col(0);
  return frame;
}

Eigen::Matrix3cd rotate_into_frame(const Eigen::Matrix3cd& alpha, const Eigen::Matrix3d& frame) {
  const Eigen::Matrix3cd R = frame.cast<std::complex<double>>();
  return R.transpose() * alpha * R;
}

std::vector<Eigen::Matrix3cd> increments(const std::vector<Eigen::Matrix3cd>& oligomer_tensors) {
  std::vector<Eigen::Matrix3cd> out;
  for (std::size_t m = 1; m < oligomer_tensors.size(); ++m)
    out.push_back(oligomer_tensors[m] - oligomer_tensors[m - 1]);
  return out;
}

UniaxialPolarizability uniaxial(const Eigen::Matrix3cd& alpha) {
  return {alpha(2, 2), 0.5 * (alpha(0, 0) + alpha(1, 1))};
}

OrientedPolarizability oriented(const UniaxialPolarizability& alpha, double order_parameter) {
  if (order_parameter < -0.5 || order_parameter > 1.0)
    throw std::invalid_argument("oriented: the order parameter must lie between -1/2 and 1");
  return {alpha.mean() - order_parameter / 3.0 * alpha.anisotropy(),
          alpha.mean() + 2.0 * order_parameter / 3.0 * alpha.anisotropy()};
}

OpticalConstants optical_constants(double wavelength, const Eigen::Matrix3cd& alpha_unit_au, double density,
                                   double order_parameter) {
  const UniaxialPolarizability u = uniaxial(alpha_unit_au);
  const OrientedPolarizability o = oriented(u, order_parameter);
  OpticalConstants result;
  result.wavelength_nm = wavelength;
  result.isotropic = lorentz_lorenz(angstrom3_per_au * u.mean(), density);
  result.ordinary = lorentz_lorenz(angstrom3_per_au * o.ordinary, density);
  result.extraordinary = lorentz_lorenz(angstrom3_per_au * o.extraordinary, density);
  return result;
}

}  // namespace ellidft::optics
