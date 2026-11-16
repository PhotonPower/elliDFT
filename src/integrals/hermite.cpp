#include "integrals/hermite.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "integrals/boys.hpp"

namespace ellidft {

HermiteE::HermiteE(int imax, int jmax, double a, double b, double A, double B)
    : imax_(imax), jmax_(jmax), tdim_(imax + jmax + 1),
      data_(static_cast<std::size_t>(imax + 1) * static_cast<std::size_t>(jmax + 1) *
                static_cast<std::size_t>(imax + jmax + 1),
            0.0) {
  const double p = a + b;
  const double AB = A - B;
  // P - A = -b/p (A - B) and P - B = a/p (A - B); written this way there is
  // no cancellation for nearby centres.
  const double PA = -b / p * AB;
  const double PB = a / p * AB;
  const double half_inverse_p = 0.5 / p;

  // Hernández Vera 2025, Eq. 19.
  data_[index(0, 0, 0)] = std::exp(-a * b / p * AB * AB);

  // Eq. 17: E_t^{i+1,j} = E_{t-1}^{ij} / (2p) + (P-A) E_t^{ij} + (t+1) E_{t+1}^{ij}
  // (the middle term is misprinted as E^t_{ij} in the source).
  for (int i = 0; i < imax; ++i)
    for (int t = 0; t <= i + 1; ++t)
      data_[index(i + 1, 0, t)] =
          half_inverse_p * (*this)(i, 0, t - 1) + PA * (*this)(i, 0, t) + (t + 1) * (*this)(i, 0, t + 1);

  // Eq. 18: E_t^{i,j+1} = E_{t-1}^{ij} / (2p) + (P-B) E_t^{ij} + (t+1) E_{t+1}^{ij}
  for (int i = 0; i <= imax; ++i)
    for (int j = 0; j < jmax; ++j)
      for (int t = 0; t <= i + j + 1; ++t)
        data_[index(i, j + 1, t)] =
            half_inverse_p * (*this)(i, j, t - 1) + PB * (*this)(i, j, t) + (t + 1) * (*this)(i, j, t + 1);
}

std::vector<HermiteIndex> hermite_indices(int L) {
  std::vector<HermiteIndex> out;
  out.reserve(static_cast<std::size_t>(n_hermite(L)));
  for (int t = 0; t <= L; ++t)
    for (int u = 0; u <= L - t; ++u)
      for (int v = 0; v <= L - t - u; ++v) out.push_back({t, u, v});
  return out;
}

HermiteCoulomb::HermiteCoulomb(int max_L) : max_L_(max_L) {
  if (max_L < 0 || max_L > boys_max_order)
    throw std::out_of_range("HermiteCoulomb: total angular momentum " + std::to_string(max_L) +
                            " is not supported");
  const auto d = static_cast<std::size_t>(max_L + 1);
  work_.assign(d * d * d * d, 0.0);
  boys_.assign(d, 0.0);
  level0_ = work_.data();
}

void HermiteCoulomb::compute(int L, double alpha, double X, double Y, double Z) {
  assert(L >= 0 && L <= max_L_);
  dim_ = L + 1;
  const auto d = static_cast<std::size_t>(dim_);
  const std::size_t cube = d * d * d;
  auto at = [&](int n, int t, int u, int v) -> double& {
    return work_[static_cast<std::size_t>(n) * cube + cube_index(t, u, v)];
  };

  // Source terms R^n_000 = (-2 alpha)^n F_n(alpha R^2) (text below Eq. 50).
  boys(L, alpha * (X * X + Y * Y + Z * Z), boys_.data());
  double factor = 1.0;
  for (int n = 0; n <= L; ++n) {
    at(n, 0, 0, 0) = factor * boys_[static_cast<std::size_t>(n)];
    factor *= -2.0 * alpha;
  }

  // Eq. 51, with the index shifted by one:
  //   R^n_{t,u,v} = (t-1) R^{n+1}_{t-2,u,v} + X R^{n+1}_{t-1,u,v}   (same for u, v).
  // Entries of total order N = t + u + v at level n need order N - 1 and
  // N - 2 at level n + 1, so the orders are built up one after the other.
  for (int N = 1; N <= L; ++N) {
    for (int n = 0; n <= L - N; ++n) {
      for (int t = 0; t <= N; ++t) {
        for (int u = 0; u <= N - t; ++u) {
          const int v = N - t - u;
          double value;
          if (t > 0) {
            value = X * at(n + 1, t - 1, u, v);
            if (t > 1) value += (t - 1) * at(n + 1, t - 2, u, v);
          } else if (u > 0) {
            value = Y * at(n + 1, t, u - 1, v);
            if (u > 1) value += (u - 1) * at(n + 1, t, u - 2, v);
          } else {
            value = Z * at(n + 1, t, u, v - 1);
            if (v > 1) value += (v - 1) * at(n + 1, t, u, v - 2);
          }
          at(n, t, u, v) = value;
        }
      }
    }
  }
  level0_ = work_.data();
}

}  // namespace ellidft
