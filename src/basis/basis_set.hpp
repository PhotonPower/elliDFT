#pragma once

#include <optional>
#include <string>
#include <vector>

#include "basis/basis_library.hpp"
#include "core/molecule.hpp"

namespace ellidft {

// Number of cartesian functions x^i y^j z^k with i+j+k = l.
constexpr int n_cartesian(int l) { return (l + 1) * (l + 2) / 2; }

// Number of real solid harmonics of degree l.
constexpr int n_spherical(int l) { return 2 * l + 1; }

// A contracted shell placed on an atom of a molecule.
struct Shell {
  int atom = 0;            // index into Molecule::atoms
  Vec3 center{};           // Bohr
  int l = 0;
  bool spherical = true;   // real solid harmonics (2l+1) or cartesian functions
  std::vector<double> exponents;
  std::vector<double> coefficients;  // as in the basis set file, not normalised
  // Coefficients d_k of the normalised radial function
  // R(r) = sum_k d_k exp(-exponents[k] r^2), see normalized_contraction().
  std::vector<double> norm_coefficients;

  std::size_t n_primitives() const { return exponents.size(); }
  int n_functions() const { return spherical ? n_spherical(l) : n_cartesian(l); }
};

// Integral of r^(2l+2) exp(-zeta r^2) from 0 to infinity,
// Gamma(l + 3/2) / (2 zeta^(l + 3/2)) (substitution u = zeta r^2 and the
// definition of the Gamma function).
double radial_gaussian_integral(int l, double zeta);

// Contraction coefficients d_k such that the radial function
//
//   R(r) = sum_k d_k exp(-exponents[k] r^2)
//
// satisfies integral_0^inf R(r)^2 r^(2l+2) dr = 1. Convention of the reference
// code: each primitive is first normalised in this sense, multiplied by the
// coefficient from the basis set file, and the contraction is normalised
// again. The angular part is handled by shell_transformation().
std::vector<double> normalized_contraction(int l, const std::vector<double>& exponents,
                                           const std::vector<double>& coefficients);

// The basis of a molecule: an ordered list of shells.
//
// Order: atoms as in the molecule; for each atom the shells as in
// BasisLibrary (by angular momentum, then file order). This is the order of
// the reference code, so matrices can be compared without permutation.
//
// All shells of a basis are either spherical or cartesian. s and p shells
// have the same number of functions in both cases.
class BasisSet {
 public:
  // Builds the basis for `molecule`. `spherical` overrides the choice declared
  // in the basis set file. Throws std::runtime_error if an element of the
  // molecule is missing in the library.
  static BasisSet build(const Molecule& molecule, const BasisLibrary& library,
                        std::optional<bool> spherical = std::nullopt);

  // Builds a basis from explicitly given shells (atom index, centre, l,
  // exponents and coefficients as in a basis set file). Function type and
  // normalised coefficients are set here; the order of the shells is kept.
  static BasisSet from_shells(std::string name, bool spherical, std::vector<Shell> shells);

  const std::string& name() const { return name_; }
  bool spherical() const { return spherical_; }

  const std::vector<Shell>& shells() const { return shells_; }
  std::size_t n_shells() const { return shells_.size(); }

  // Total number of basis functions.
  int n_functions() const { return n_functions_; }

  // Index of the first basis function of shell `i`.
  int first_function(std::size_t i) const { return first_function_.at(i); }

  int max_l() const { return max_l_; }
  std::size_t n_primitives() const;

 private:
  std::string name_;
  bool spherical_ = true;
  std::vector<Shell> shells_;
  std::vector<int> first_function_;
  int n_functions_ = 0;
  int max_l_ = 0;
};

}  // namespace ellidft
