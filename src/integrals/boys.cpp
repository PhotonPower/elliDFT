#include "integrals/boys.hpp"

#include <cassert>
#include <cmath>
#include <vector>

namespace ellidft {

namespace {

// --- Method ---------------------------------------------------------------
//
// x < x_switch: Taylor expansion around the nearest point x0 of an equidistant
//   table. From the definition, dF_n/dx = -F_(n+1), hence
//       F_n(x0 + d) = sum_k F_(n+k)(x0) (-d)^k / k!.
//   With |d| <= spacing/2 = 0.05 and F_(n+k) <= F_n, the remainder after
//   k = 0..7 is below 0.05^8 / 8! = 1e-15 relative. Only the highest order is
//   expanded; lower orders follow from the downward recursion
//       F_k = (2x F_(k+1) + exp(-x)) / (2k+1)
//   (Vikhamar-Sandberg & Repisky 2025, Eq. 15), which adds positive terms and
//   is stable for all x.
//
// x >= x_switch: F_0 from the asymptotic formula F_0 = sqrt(pi/x)/2
//   (Vikhamar-Sandberg & Repisky 2025, Eq. 11 with Gamma(1/2) = sqrt(pi)); its
//   relative error is erfc(sqrt(x)) < 3e-17 for x >= 36. Higher orders follow
//   from the upward recursion
//       F_(k+1) = ((2k+1) F_k - exp(-x)) / (2x)
//   (ibid., Eq. 14). The subtraction is harmless here because
//   exp(-x) << (2k+1) F_k for x >= 36 and k <= boys_max_order.

constexpr double x_switch = 36.0;
constexpr double spacing = 0.1;
constexpr int taylor_terms = 8;                                   // k = 0..7
constexpr int table_orders = boys_max_order + taylor_terms;       // orders 0..table_orders-1
constexpr int table_points = 361;                                 // x0 = 0, 0.1, ..., 36.0

struct Table {
  std::vector<double> values;  // values[point * table_orders + order]
  Table() : values(static_cast<std::size_t>(table_points) * table_orders) {
    for (int i = 0; i < table_points; ++i)
      boys_series(table_orders - 1, spacing * i,
                  &values[static_cast<std::size_t>(i) * table_orders]);
  }
};

const Table& table() {
  static const Table instance;
  return instance;
}

constexpr double inverse_factorial[taylor_terms] = {1.0,        1.0,         1.0 / 2.0,   1.0 / 6.0,
                                                    1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0, 1.0 / 5040.0};

}  // namespace

void boys_series(int nmax, double x, double* out) {
  assert(nmax >= 0 && x >= 0.0);
  // F_n(x) = exp(-x) sum_{k>=0} (2x)^k / [(2n+1)(2n+3)...(2n+2k+1)],
  // obtained by applying the downward recursion (Eq. 15 of Vikhamar-Sandberg
  // & Repisky 2025) repeatedly to F_n. All terms are positive.
  const double exp_minus_x = std::exp(-x);
  double term = 1.0 / (2 * nmax + 1);
  double sum = term;
  for (int k = 1;; ++k) {
    term *= 2.0 * x / (2 * nmax + 2 * k + 1);
    sum += term;
    // Past the maximum of the terms they decrease faster than geometrically.
    if (term < 1e-17 * sum && 2 * nmax + 2 * k + 1 > 2.0 * x) break;
  }
  out[nmax] = exp_minus_x * sum;
  for (int k = nmax - 1; k >= 0; --k)
    out[k] = (2.0 * x * out[k + 1] + exp_minus_x) / (2 * k + 1);
}

void boys(int nmax, double x, double* out) {
  assert(nmax >= 0 && nmax <= boys_max_order && x >= 0.0);

  if (x < x_switch) {
    const int point = static_cast<int>(x / spacing + 0.5);
    const double d = x - spacing * point;
    const double* f = &table().values[static_cast<std::size_t>(point) * table_orders + nmax];
    double sum = 0.0;
    double power = 1.0;  // (-d)^k
    for (int k = 0; k < taylor_terms; ++k) {
      sum += f[k] * power * inverse_factorial[k];
      power *= -d;
    }
    out[nmax] = sum;
    if (nmax > 0) {
      const double exp_minus_x = std::exp(-x);
      for (int k = nmax - 1; k >= 0; --k)
        out[k] = (2.0 * x * out[k + 1] + exp_minus_x) / (2 * k + 1);
    }
    return;
  }

  out[0] = 0.5 * std::sqrt(M_PI / x);
  if (nmax > 0) {
    const double exp_minus_x = std::exp(-x);
    const double inverse_2x = 0.5 / x;
    for (int k = 0; k < nmax; ++k) out[k + 1] = ((2 * k + 1) * out[k] - exp_minus_x) * inverse_2x;
  }
}

}  // namespace ellidft
