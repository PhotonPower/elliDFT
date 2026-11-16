#pragma once

// Explicit matrices A and B of linear response theory, built term by term
// from the full integral tensor. For tests with small molecules only; the
// library applies (A+B) and (A-B) without building them (specification,
// section 4).

#include <vector>

#include <Eigen/Dense>

#include "full_eri.hpp"
#include "scf/rhf.hpp"
#include "scf/xc_matrix.hpp"

namespace ellidft::testing {

struct ExplicitAB {
  Eigen::MatrixXd A, B;
};

// (pq|rs) over molecular orbitals from the AO tensor, index ((p*n+q)*n+r)*n+s.
inline FullEri transform_to_mo(const FullEri& ao, const Eigen::MatrixXd& C) {
  const int n = ao.n;
  const int m = static_cast<int>(C.cols());
  // One index at a time.
  auto quarter = [&](const std::vector<double>& in, int d0, int d1, int d2, int d3, int position) {
    const int dims[4] = {d0, d1, d2, d3};
    int out_dims[4] = {d0, d1, d2, d3};
    out_dims[position] = m;
    std::vector<double> out(static_cast<std::size_t>(out_dims[0]) * out_dims[1] * out_dims[2] * out_dims[3], 0.0);
    int idx[4];
    for (idx[0] = 0; idx[0] < dims[0]; ++idx[0])
      for (idx[1] = 0; idx[1] < dims[1]; ++idx[1])
        for (idx[2] = 0; idx[2] < dims[2]; ++idx[2])
          for (idx[3] = 0; idx[3] < dims[3]; ++idx[3]) {
            const double value = in[((static_cast<std::size_t>(idx[0]) * dims[1] + idx[1]) * dims[2] + idx[2]) * dims[3] + idx[3]];
            if (value == 0.0) continue;
            int o[4] = {idx[0], idx[1], idx[2], idx[3]};
            for (int p = 0; p < m; ++p) {
              o[position] = p;
              out[((static_cast<std::size_t>(o[0]) * out_dims[1] + o[1]) * out_dims[2] + o[2]) * out_dims[3] + o[3]] +=
                  C(idx[position], p) * value;
            }
          }
    return out;
  };
  std::vector<double> t = quarter(ao.values, n, n, n, n, 0);
  t = quarter(t, m, n, n, n, 1);
  t = quarter(t, m, m, n, n, 2);
  t = quarter(t, m, m, m, n, 3);
  FullEri mo;
  mo.n = m;
  mo.values = std::move(t);
  return mo;
}

// A_{ia,jb} = delta (eps_a - eps_i) + 2 (ia|jb) + 2 (ia|f_xc|jb) - c_x (ij|ab) - c_lr (ij|ab)_erf
// B_{ia,jb} =                         2 (ia|jb) + 2 (ia|f_xc|jb) - c_x (ib|ja) - c_lr (ib|ja)_erf
// (specification, section 8). Pair index i * n_virt + a.
inline ExplicitAB explicit_ab(const Molecule& molecule, const BasisSet& basis, const Method& method,
                              const ScfResult& scf) {
  const Eigen::MatrixXd& C = scf.coefficients;
  const int n_occ = scf.n_occupied;
  const int n_virt = static_cast<int>(C.cols()) - n_occ;
  const int N = n_occ * n_virt;
  const Functional& f = method.functional;

  const FullEri coulomb = transform_to_mo(full_eri_tensor(basis), C);
  FullEri attenuated;
  if (f.long_range_exchange != 0.0) attenuated = transform_to_mo(full_eri_tensor(basis, f.omega), C);

  // (ia|f_xc|jb): kernel applied to the orbital products of all pairs (jb).
  std::vector<Eigen::MatrixXd> kernel;
  if (f.has_density_part()) {
    const auto xc = make_libxc_evaluator(f);
    const MolecularGrid grid = build_molecular_grid(molecule, method.grid);
    std::vector<Eigen::MatrixXd> products;
    for (int j = 0; j < n_occ; ++j)
      for (int b = 0; b < n_virt; ++b) products.push_back(C.col(j) * C.col(n_occ + b).transpose());
    kernel = integrate_xc_kernel(basis, grid, *xc, scf.density, products);
  }

  ExplicitAB result;
  result.A = Eigen::MatrixXd::Zero(N, N);
  result.B = Eigen::MatrixXd::Zero(N, N);
  for (int i = 0; i < n_occ; ++i)
    for (int a = 0; a < n_virt; ++a)
      for (int j = 0; j < n_occ; ++j)
        for (int b = 0; b < n_virt; ++b) {
          const int ia = i * n_virt + a, jb = j * n_virt + b;
          const int A_ = n_occ + a, B_ = n_occ + b;
          double coupling = 2.0 * coulomb(i, A_, j, B_);
          if (!kernel.empty())
            coupling += 2.0 * C.col(i).dot(kernel[static_cast<std::size_t>(jb)] * C.col(A_));
          double a_value = coupling - f.exact_exchange * coulomb(i, j, A_, B_);
          double b_value = coupling - f.exact_exchange * coulomb(i, B_, j, A_);
          if (f.long_range_exchange != 0.0) {
            a_value -= f.long_range_exchange * attenuated(i, j, A_, B_);
            b_value -= f.long_range_exchange * attenuated(i, B_, j, A_);
          }
          if (ia == jb) a_value += scf.orbital_energies(A_) - scf.orbital_energies(i);
          result.A(ia, jb) = a_value;
          result.B(ia, jb) = b_value;
        }
  return result;
}

}  // namespace ellidft::testing
