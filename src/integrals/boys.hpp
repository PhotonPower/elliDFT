#pragma once

namespace ellidft {

// Highest order n for which boys() can be called. Four-centre integrals over
// shells with l <= 6 need n <= 24; the rest is margin.
inline constexpr int boys_max_order = 32;

// Boys function
//
//   F_n(x) = integral_0^1 t^(2n) exp(-x t^2) dt,   x >= 0
//
// (Vikhamar-Sandberg & Repisky 2025, Eq. 3). Writes F_0(x) .. F_nmax(x) to
// out[0..nmax]. Relative accuracy about 1e-14 for all n <= boys_max_order.
//
// Requires 0 <= nmax <= boys_max_order and x >= 0 (not checked in release
// builds; this is the innermost function of the integral code).
void boys(int nmax, double x, double* out);

// Reference implementation used to fill the interpolation table and in tests:
// series with positive terms for F_nmax, then downward recursion. Accurate for
// all x >= 0 but slow for large x (the number of terms grows with x).
void boys_series(int nmax, double x, double* out);

}  // namespace ellidft
