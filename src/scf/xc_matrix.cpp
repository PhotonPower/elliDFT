#include "scf/xc_matrix.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "core/parallel.hpp"
#include "grid/basis_on_grid.hpp"

namespace ellidft {

XcIntegrals integrate_xc(const BasisSet& basis, const MolecularGrid& grid, const XcEvaluator& functional,
                         const Eigen::MatrixXd& density, std::size_t block_size, double basis_cutoff) {
  const Eigen::Index n = basis.n_functions();
  if (density.rows() != n || density.cols() != n)
    throw std::invalid_argument("integrate_xc: density matrix has the wrong dimension");
  if (block_size == 0) throw std::invalid_argument("integrate_xc: block size must be positive");

  const bool gga = functional.needs_gradient();
  const Eigen::MatrixXd D = 0.5 * (density + density.transpose());

  XcIntegrals result;
  result.matrix = Eigen::MatrixXd::Zero(n, n);

  // Blocks are distributed over the threads; every thread accumulates its own
  // sums, which are added in thread order afterwards.
  const long n_blocks = static_cast<long>((grid.size() + block_size - 1) / block_size);
  const int n_threads = max_threads();
  std::vector<Eigen::MatrixXd> half_t(static_cast<std::size_t>(n_threads));  // V = half + half^T
  std::vector<double> energy_t(static_cast<std::size_t>(n_threads), 0.0), electrons_t(energy_t.size(), 0.0);

#pragma omp parallel num_threads(n_threads)
  {
    const auto thread = static_cast<std::size_t>(thread_index());
    Eigen::MatrixXd& half = half_t[thread];
    half = Eigen::MatrixXd::Zero(n, n);
    double energy = 0.0, electrons = 0.0;

    BasisValues phi;
    std::vector<int> active;  // basis functions that are not negligible in the block
    Eigen::MatrixXd D_active, phiD, B, half_active;
    std::vector<double> rho, sigma, gx, gy, gz, exc, vrho, vsigma;

#pragma omp for schedule(dynamic, 4)
    for (long block = 0; block < n_blocks; ++block) {
      const std::size_t start = static_cast<std::size_t>(block) * block_size;
      const std::size_t count = std::min(block_size, grid.size() - start);
      const auto np = static_cast<Eigen::Index>(count);
      evaluate_basis_active(basis, &grid.points[start], count, gga, phi, active, basis_cutoff);
      const auto na = static_cast<Eigen::Index>(active.size());
      if (na == 0) continue;
      D_active.resize(na, na);
      for (Eigen::Index i = 0; i < na; ++i)
        for (Eigen::Index j = 0; j < na; ++j)
          D_active(i, j) = D(active[static_cast<std::size_t>(i)], active[static_cast<std::size_t>(j)]);

      // rho_p = sum_mu,nu chi_mu(p) D_mu,nu chi_nu(p);
      // grad rho_p = 2 sum_mu,nu grad chi_mu(p) D_mu,nu chi_nu(p) for symmetric D.
      phiD.noalias() = phi.value * D_active;
      rho.resize(count);
      for (Eigen::Index p = 0; p < np; ++p) rho[static_cast<std::size_t>(p)] = phiD.row(p).dot(phi.value.row(p));
      if (gga) {
        sigma.resize(count);
        gx.resize(count);
        gy.resize(count);
        gz.resize(count);
        for (Eigen::Index p = 0; p < np; ++p) {
          const auto i = static_cast<std::size_t>(p);
          gx[i] = 2.0 * phiD.row(p).dot(phi.dx.row(p));
          gy[i] = 2.0 * phiD.row(p).dot(phi.dy.row(p));
          gz[i] = 2.0 * phiD.row(p).dot(phi.dz.row(p));
          sigma[i] = gx[i] * gx[i] + gy[i] * gy[i] + gz[i] * gz[i];
        }
      }

      exc.resize(count);
      vrho.resize(count);
      vsigma.resize(count);
      functional.evaluate(count, rho.data(), gga ? sigma.data() : nullptr, exc.data(), vrho.data(), vsigma.data());

      // Chain rule: dE = integral (vrho d rho + vsigma d sigma), with
      // d rho = chi_mu chi_nu dD and d sigma = 2 grad rho . grad(chi_mu chi_nu) dD:
      //   V_mu,nu = sum_p w_p [ vrho chi_mu chi_nu
      //                         + 2 vsigma grad rho . (grad chi_mu chi_nu + chi_mu grad chi_nu) ].
      // With B(p, mu) = w_p (vrho chi_mu / 2 + 2 vsigma grad rho . grad chi_mu)
      // this is V = chi^T B + B^T chi.
      B.resize(np, na);
      for (Eigen::Index p = 0; p < np; ++p) {
        const auto i = static_cast<std::size_t>(p);
        const double w = grid.weights[start + i];
        energy += w * exc[i];
        electrons += w * rho[i];
        B.row(p) = (0.5 * w * vrho[i]) * phi.value.row(p);
        if (gga) {
          const double f = 2.0 * w * vsigma[i];
          B.row(p) += f * (gx[i] * phi.dx.row(p) + gy[i] * phi.dy.row(p) + gz[i] * phi.dz.row(p));
        }
      }
      half_active.noalias() = phi.value.transpose() * B;
      for (Eigen::Index i = 0; i < na; ++i)
        for (Eigen::Index j = 0; j < na; ++j)
          half(active[static_cast<std::size_t>(i)], active[static_cast<std::size_t>(j)]) += half_active(i, j);
    }
    energy_t[thread] = energy;
    electrons_t[thread] = electrons;
  }

  Eigen::MatrixXd half = Eigen::MatrixXd::Zero(n, n);
  for (std::size_t t = 0; t < half_t.size(); ++t) {
    if (half_t[t].size() == 0) continue;
    half += half_t[t];
    result.energy += energy_t[t];
    result.n_electrons += electrons_t[t];
  }
  result.matrix = half + half.transpose();
  return result;
}

std::vector<Eigen::MatrixXd> integrate_xc_kernel(const BasisSet& basis, const MolecularGrid& grid,
                                                 const XcEvaluator& functional, const Eigen::MatrixXd& density,
                                                 const std::vector<Eigen::MatrixXd>& trial_densities,
                                                 std::size_t block_size, double basis_cutoff) {
  const Eigen::Index n = basis.n_functions();
  if (density.rows() != n || density.cols() != n)
    throw std::invalid_argument("integrate_xc_kernel: density matrix has the wrong dimension");
  for (const Eigen::MatrixXd& P : trial_densities)
    if (P.rows() != n || P.cols() != n)
      throw std::invalid_argument("integrate_xc_kernel: trial density has the wrong dimension");
  if (block_size == 0) throw std::invalid_argument("integrate_xc_kernel: block size must be positive");

  const bool gga = functional.needs_gradient();
  const std::size_t n_trial = trial_densities.size();
  const Eigen::MatrixXd D = 0.5 * (density + density.transpose());
  std::vector<Eigen::MatrixXd> P_sym;
  for (const Eigen::MatrixXd& P : trial_densities) P_sym.push_back(0.5 * (P + P.transpose()));

  const long n_blocks = static_cast<long>((grid.size() + block_size - 1) / block_size);
  const int n_threads = max_threads();
  const Eigen::MatrixXd zero = Eigen::MatrixXd::Zero(n, n);
  std::vector<std::vector<Eigen::MatrixXd>> half_t(static_cast<std::size_t>(n_threads));

#pragma omp parallel num_threads(n_threads)
  {
    std::vector<Eigen::MatrixXd>& half = half_t[static_cast<std::size_t>(thread_index())];
    half.assign(n_trial, zero);

    BasisValues phi;
    std::vector<int> active;
    Eigen::MatrixXd D_active, P_active, phiD, phiP, B, half_active;
    std::vector<double> rho, sigma, gx, gy, gz, vrho, vsigma, v2rho2, v2rhosigma, v2sigma2;

#pragma omp for schedule(dynamic, 4)
    for (long block = 0; block < n_blocks; ++block) {
      const std::size_t start = static_cast<std::size_t>(block) * block_size;
      const std::size_t count = std::min(block_size, grid.size() - start);
      const auto np = static_cast<Eigen::Index>(count);
      evaluate_basis_active(basis, &grid.points[start], count, gga, phi, active, basis_cutoff);
      const auto na = static_cast<Eigen::Index>(active.size());
      if (na == 0) continue;
      auto restrict_to_active = [&](const Eigen::MatrixXd& full, Eigen::MatrixXd& part) {
        part.resize(na, na);
        for (Eigen::Index i = 0; i < na; ++i)
          for (Eigen::Index j = 0; j < na; ++j)
            part(i, j) = full(active[static_cast<std::size_t>(i)], active[static_cast<std::size_t>(j)]);
      };

      // Ground-state density and its gradient.
      restrict_to_active(D, D_active);
      phiD.noalias() = phi.value * D_active;
      rho.resize(count);
      sigma.assign(count, 0.0);
      gx.assign(count, 0.0);
      gy.assign(count, 0.0);
      gz.assign(count, 0.0);
      for (Eigen::Index p = 0; p < np; ++p) {
        const auto i = static_cast<std::size_t>(p);
        rho[i] = phiD.row(p).dot(phi.value.row(p));
        if (gga) {
          gx[i] = 2.0 * phiD.row(p).dot(phi.dx.row(p));
          gy[i] = 2.0 * phiD.row(p).dot(phi.dy.row(p));
          gz[i] = 2.0 * phiD.row(p).dot(phi.dz.row(p));
          sigma[i] = gx[i] * gx[i] + gy[i] * gy[i] + gz[i] * gz[i];
        }
      }
      vrho.resize(count);
      vsigma.resize(count);
      v2rho2.resize(count);
      v2rhosigma.resize(count);
      v2sigma2.resize(count);
      functional.evaluate_kernel(count, rho.data(), gga ? sigma.data() : nullptr, vrho.data(), vsigma.data(),
                                 v2rho2.data(), v2rhosigma.data(), v2sigma2.data());

      // Variation of V_mu,nu = integral e_rho chi_mu chi_nu + 2 e_sigma grad rho . grad(chi_mu chi_nu)
      // for a first-order density rho1 with sigma1 = 2 grad rho . grad rho1:
      //   dV_mu,nu = integral (e_rr rho1 + e_rs sigma1) chi_mu chi_nu
      //            + [2 (e_rs rho1 + e_ss sigma1) grad rho + 2 e_s grad rho1] . grad(chi_mu chi_nu),
      // where e_r, e_s, e_rr, e_rs, e_ss are the derivatives of the energy per
      // volume with respect to rho and sigma. With
      //   B(p, mu) = w_p [ a chi_mu / 2 + g . grad chi_mu ],  a = e_rr rho1 + e_rs sigma1,
      //   g = 2 (e_rs rho1 + e_ss sigma1) grad rho + 2 e_s grad rho1,
      // this is dV = chi^T B + B^T chi.
      for (std::size_t t = 0; t < n_trial; ++t) {
        restrict_to_active(P_sym[t], P_active);
        phiP.noalias() = phi.value * P_active;
        B.resize(np, na);
        for (Eigen::Index p = 0; p < np; ++p) {
          const auto i = static_cast<std::size_t>(p);
          const double w = grid.weights[start + i];
          const double rho1 = phiP.row(p).dot(phi.value.row(p));
          if (!gga) {
            B.row(p) = (0.5 * w * v2rho2[i] * rho1) * phi.value.row(p);
            continue;
          }
          const double g1x = 2.0 * phiP.row(p).dot(phi.dx.row(p));
          const double g1y = 2.0 * phiP.row(p).dot(phi.dy.row(p));
          const double g1z = 2.0 * phiP.row(p).dot(phi.dz.row(p));
          const double sigma1 = 2.0 * (gx[i] * g1x + gy[i] * g1y + gz[i] * g1z);
          const double a = v2rho2[i] * rho1 + v2rhosigma[i] * sigma1;
          const double b = 2.0 * (v2rhosigma[i] * rho1 + v2sigma2[i] * sigma1);
          const double c = 2.0 * vsigma[i];
          B.row(p) = (0.5 * w * a) * phi.value.row(p) + (w * (b * gx[i] + c * g1x)) * phi.dx.row(p) +
                     (w * (b * gy[i] + c * g1y)) * phi.dy.row(p) + (w * (b * gz[i] + c * g1z)) * phi.dz.row(p);
        }
        half_active.noalias() = phi.value.transpose() * B;
        for (Eigen::Index i = 0; i < na; ++i)
          for (Eigen::Index j = 0; j < na; ++j)
            half[t](active[static_cast<std::size_t>(i)], active[static_cast<std::size_t>(j)]) += half_active(i, j);
      }
    }
  }

  std::vector<Eigen::MatrixXd> result(n_trial, zero);
  for (std::size_t t = 0; t < n_trial; ++t) {
    Eigen::MatrixXd half = zero;
    for (const auto& thread_half : half_t)
      if (!thread_half.empty()) half += thread_half[t];
    result[t] = half + half.transpose();
  }
  return result;
}

}  // namespace ellidft
