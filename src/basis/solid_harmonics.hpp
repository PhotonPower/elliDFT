#pragma once

#include <vector>

#include <Eigen/Dense>

namespace ellidft {

// Exponents of a cartesian monomial x^x y^y z^z.
struct CartesianPowers {
  int x = 0, y = 0, z = 0;
};

// The (l+1)(l+2)/2 cartesian monomials of degree l in the order of the
// reference code (lexicographic): xx, xy, xz, yy, yz, zz for l = 2.
const std::vector<CartesianPowers>& cartesian_powers(int l);

// Position of x^lx y^ly z^lz in cartesian_powers(lx + ly + lz).
constexpr int cartesian_index([[maybe_unused]] int lx, int ly, int lz) {
  const int rest = ly + lz;  // components with larger lx come first
  return rest * (rest + 1) / 2 + lz;
}

// Order of the real solid harmonics in a spherical shell, as values of m.
// Reference code convention: p as x, y, z, i.e. m = +1, -1, 0; all other l as
// m = -l .. +l.
std::vector<int> spherical_m_order(int l);

// Unnormalised real solid harmonic of degree l as a polynomial: coefficients
// of the monomials in cartesian_powers(l) order.
//
//   m >= 0:  Re[ r^l P_l^|m|(cos theta) exp(i |m| phi) ]
//   m <  0:  Im[ r^l P_l^|m|(cos theta) exp(i |m| phi) ]
//
// with the associated Legendre function without Condon-Shortley phase
// (Ribaldone & Desmarais 2025, Eqs. 1, 2, 5, 14).
std::vector<double> real_solid_harmonic(int l, int m);

// Integral of x^a y^b z^c over the unit sphere.
double sphere_monomial_integral(int a, int b, int c);

// Matrix that turns the cartesian functions of a shell into its basis
// functions: chi_i = sum_c T(i, c) * x^cx y^cy z^cz * R(r), where R(r) is the
// contracted radial function with integral R^2 r^(2l+2) dr = 1 (see
// normalized_contraction in basis_set.hpp).
//
//  * spherical: rows are the real solid harmonics in spherical_m_order(l),
//    normalised so that every basis function has norm 1.
//  * cartesian: for l <= 1 the same functions as in the spherical case
//    (s; px, py, pz). For l >= 2 the identity: the functions carry only the
//    radial normalisation, as in the reference code (libcint). Their norm is
//    not 1, e.g. <xx|xx> = 4 pi / 5.
const Eigen::MatrixXd& shell_transformation(int l, bool spherical);

}  // namespace ellidft
