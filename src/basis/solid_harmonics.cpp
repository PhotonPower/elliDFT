#include "basis/solid_harmonics.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include "basis/basis_library.hpp"

namespace ellidft {

namespace {

void check_l(int l) {
  if (l < 0 || l > max_angular_momentum)
    throw std::out_of_range("angular momentum l = " + std::to_string(l) + " is not supported");
}

double factorial(int n) {
  double f = 1.0;
  for (int k = 2; k <= n; ++k) f *= k;
  return f;
}

double binomial(int n, int k) { return factorial(n) / (factorial(k) * factorial(n - k)); }

std::array<std::vector<CartesianPowers>, max_angular_momentum + 1> make_cartesian_powers() {
  std::array<std::vector<CartesianPowers>, max_angular_momentum + 1> out;
  for (int l = 0; l <= max_angular_momentum; ++l)
    for (int lx = l; lx >= 0; --lx)
      for (int ly = l - lx; ly >= 0; --ly) out[static_cast<std::size_t>(l)].push_back({lx, ly, l - lx - ly});
  return out;
}

Eigen::MatrixXd make_spherical_transformation(int l) {
  const auto& powers = cartesian_powers(l);
  const auto n_cart = static_cast<Eigen::Index>(powers.size());
  const std::vector<int> m_order = spherical_m_order(l);
  Eigen::MatrixXd T(static_cast<Eigen::Index>(m_order.size()), n_cart);

  for (Eigen::Index row = 0; row < T.rows(); ++row) {
    const std::vector<double> polynomial = real_solid_harmonic(l, m_order[static_cast<std::size_t>(row)]);
    // Norm on the unit sphere: integral of (sum_c a_c x^.. y^.. z^..)^2.
    double norm2 = 0.0;
    for (Eigen::Index c = 0; c < n_cart; ++c)
      for (Eigen::Index d = 0; d < n_cart; ++d) {
        const auto& pc = powers[static_cast<std::size_t>(c)];
        const auto& pd = powers[static_cast<std::size_t>(d)];
        norm2 += polynomial[static_cast<std::size_t>(c)] * polynomial[static_cast<std::size_t>(d)] *
                 sphere_monomial_integral(pc.x + pd.x, pc.y + pd.y, pc.z + pd.z);
      }
    for (Eigen::Index c = 0; c < n_cart; ++c)
      T(row, c) = polynomial[static_cast<std::size_t>(c)] / std::sqrt(norm2);
  }
  return T;
}

struct Transformations {
  std::array<Eigen::MatrixXd, max_angular_momentum + 1> spherical, cartesian;
  Transformations() {
    for (int l = 0; l <= max_angular_momentum; ++l) {
      const auto i = static_cast<std::size_t>(l);
      spherical[i] = make_spherical_transformation(l);
      if (l <= 1) {
        cartesian[i] = spherical[i];
      } else {
        const auto n = static_cast<Eigen::Index>(cartesian_powers(l).size());
        cartesian[i] = Eigen::MatrixXd::Identity(n, n);
      }
    }
  }
};

}  // namespace

const std::vector<CartesianPowers>& cartesian_powers(int l) {
  check_l(l);
  static const auto table = make_cartesian_powers();
  return table[static_cast<std::size_t>(l)];
}

std::vector<int> spherical_m_order(int l) {
  check_l(l);
  if (l == 1) return {1, -1, 0};
  std::vector<int> order;
  for (int m = -l; m <= l; ++m) order.push_back(m);
  return order;
}

std::vector<double> real_solid_harmonic(int l, int m) {
  check_l(l);
  const int M = m < 0 ? -m : m;
  if (M > l) throw std::out_of_range("real_solid_harmonic: |m| > l");

  std::vector<double> coefficients(cartesian_powers(l).size(), 0.0);

  // Ribaldone & Desmarais 2025, Eq. 14:
  //   Y_l^M = 1/(2^l l!) (x+iy)^M sum_{k=0}^{floor((l-M)/2)} binom(l,k) (-1)^k
  //           (2l-2k)! / (l-M-2k)!  (x^2+y^2+z^2)^k  z^(l-M-2k)
  // The real part (m >= 0) or imaginary part (m < 0) is taken, Eq. 5.
  const double prefactor = 1.0 / (std::pow(2.0, l) * factorial(l));
  for (int k = 0; 2 * k <= l - M; ++k) {
    const double a_k = prefactor * binomial(l, k) * ((k % 2 == 0) ? 1.0 : -1.0) *
                       factorial(2 * l - 2 * k) / factorial(l - M - 2 * k);
    // (x^2+y^2+z^2)^k = sum_{a+b+c=k} k!/(a! b! c!) x^2a y^2b z^2c
    for (int a = 0; a <= k; ++a) {
      for (int b = 0; b <= k - a; ++b) {
        const int c = k - a - b;
        const double multinomial = factorial(k) / (factorial(a) * factorial(b) * factorial(c));
        // (x+iy)^M = sum_s binom(M,s) x^(M-s) (iy)^s; i^s is real for even s
        // (sign (-1)^(s/2)) and imaginary for odd s (sign (-1)^((s-1)/2)).
        for (int s = 0; s <= M; ++s) {
          const bool real_term = (s % 2 == 0);
          if (real_term != (m >= 0)) continue;
          const double sign = ((s / 2) % 2 == 0) ? 1.0 : -1.0;
          const int px = 2 * a + M - s;
          const int py = 2 * b + s;
          const int pz = 2 * c + l - M - 2 * k;
          coefficients[static_cast<std::size_t>(cartesian_index(px, py, pz))] +=
              a_k * multinomial * binomial(M, s) * sign;
        }
      }
    }
  }
  return coefficients;
}

double sphere_monomial_integral(int a, int b, int c) {
  if (a < 0 || b < 0 || c < 0) throw std::out_of_range("sphere_monomial_integral: negative power");
  if (a % 2 != 0 || b % 2 != 0 || c % 2 != 0) return 0.0;
  // Integrate x^a y^b z^c exp(-r^2) over all space in two ways. Cartesian:
  // product of three integrals, each Gamma((n+1)/2) for even n. Spherical:
  // radial integral Gamma((a+b+c+3)/2) / 2 times the integral over the unit
  // sphere. Both follow from the definition of the Gamma function.
  return 2.0 * std::tgamma(0.5 * (a + 1)) * std::tgamma(0.5 * (b + 1)) * std::tgamma(0.5 * (c + 1)) /
         std::tgamma(0.5 * (a + b + c + 3));
}

const Eigen::MatrixXd& shell_transformation(int l, bool spherical) {
  check_l(l);
  static const Transformations table;
  return spherical ? table.spherical[static_cast<std::size_t>(l)]
                   : table.cartesian[static_cast<std::size_t>(l)];
}

}  // namespace ellidft
