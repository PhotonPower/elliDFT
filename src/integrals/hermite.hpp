#pragma once

// Building blocks of the McMurchie-Davidson scheme: Hermite expansion
// coefficients of products of cartesian Gaussians and Hermite Coulomb
// integrals. Equations from Hernández Vera 2025 (arXiv:2510.24968), see
// docs/sources/README.md.

#include <cassert>
#include <vector>

namespace ellidft {

// Hermite expansion coefficients in one cartesian direction for one pair of
// primitives (Hernández Vera 2025, Eq. 16, with the factor exp(-mu X_AB^2)
// included as in Eq. 33):
//
//   (x-A)^i (x-B)^j exp(-a (x-A)^2) exp(-b (x-B)^2) = sum_{t=0}^{i+j} E_t^{ij} Lambda_t(x-P)
//
// with p = a + b, P = (a A + b B) / p and the Hermite Gaussian
// Lambda_t(x-P) = (d/dP)^t exp(-p (x-P)^2) (Eq. 13).
class HermiteE {
 public:
  // Coefficients for 0 <= i <= imax, 0 <= j <= jmax, 0 <= t <= i + j.
  HermiteE(int imax, int jmax, double a, double b, double A, double B);

  // E_t^{ij}; zero for t < 0 or t > i + j (Eq. 20).
  double operator()(int i, int j, int t) const {
    assert(i >= 0 && i <= imax_ && j >= 0 && j <= jmax_);
    if (t < 0 || t > i + j) return 0.0;
    return data_[index(i, j, t)];
  }

 private:
  std::size_t index(int i, int j, int t) const {
    return (static_cast<std::size_t>(i) * static_cast<std::size_t>(jmax_ + 1) +
            static_cast<std::size_t>(j)) * static_cast<std::size_t>(tdim_) +
           static_cast<std::size_t>(t);
  }
  [[maybe_unused]] int imax_;  // only read by assertions
  int jmax_, tdim_;
  std::vector<double> data_;
};

// Number of Hermite Gaussians Lambda_tuv with t + u + v <= L.
constexpr int n_hermite(int L) { return (L + 1) * (L + 2) * (L + 3) / 6; }

struct HermiteIndex {
  int t, u, v;
};

// All (t, u, v) with t + u + v <= L, in a fixed order (by t, then u, then v).
std::vector<HermiteIndex> hermite_indices(int L);

// Hermite Coulomb integrals R^0_{tuv}(alpha, R) for all t + u + v <= L
// (Hernández Vera 2025, Eq. 50), from the Boys function and the recurrences
// of Eq. 51. Reusable workspace; not thread safe, use one object per thread.
class HermiteCoulomb {
 public:
  explicit HermiteCoulomb(int max_L);

  // Computes the integrals for the given L <= max_L. Afterwards
  // value(t, u, v) is R^0_{tuv} and data()[cube_index(t, u, v)] likewise,
  // where the cube has edge length L + 1.
  void compute(int L, double alpha, double X, double Y, double Z);

  double value(int t, int u, int v) const { return level0_[cube_index(t, u, v)]; }
  const double* data() const { return level0_; }
  std::size_t cube_index(int t, int u, int v) const {
    const auto d = static_cast<std::size_t>(dim_);
    return (static_cast<std::size_t>(t) * d + static_cast<std::size_t>(u)) * d +
           static_cast<std::size_t>(v);
  }
  int dim() const { return dim_; }

 private:
  [[maybe_unused]] int max_L_;  // only read by assertions
  int dim_ = 1;                // L + 1 of the last call
  std::vector<double> work_;   // work_[n * dim^3 + cube index] = R^n_{tuv}
  std::vector<double> boys_;
  const double* level0_ = nullptr;
};

}  // namespace ellidft
