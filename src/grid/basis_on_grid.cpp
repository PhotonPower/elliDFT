#include "grid/basis_on_grid.hpp"

#include <algorithm>
#include <cmath>

#include "basis/solid_harmonics.hpp"

namespace ellidft {

namespace {

// Upper bound of |sum_k d_k r^l exp(-a_k r^2)| for all r >= r_min. Each term
// r^l exp(-a r^2) has its maximum at r = sqrt(l / (2a)) and decreases beyond.
double radial_bound(const Shell& shell, double r_min) {
  double bound = 0.0;
  for (std::size_t k = 0; k < shell.n_primitives(); ++k) {
    const double a = shell.exponents[k];
    const double r = std::max(r_min, std::sqrt(shell.l / (2.0 * a)));
    bound += std::abs(shell.norm_coefficients[k]) * std::pow(r, shell.l) * std::exp(-a * r * r);
  }
  return bound;
}

}  // namespace

void evaluate_basis(const BasisSet& basis, const Vec3* points, std::size_t n_points, bool gradients,
                    BasisValues& out, double cutoff) {
  BasisValues compact;
  std::vector<int> active;
  evaluate_basis_active(basis, points, n_points, gradients, compact, active, cutoff);
  const auto np = static_cast<Eigen::Index>(n_points);
  const Eigen::Index nf = basis.n_functions();
  out.value.setZero(np, nf);
  if (gradients) {
    out.dx.setZero(np, nf);
    out.dy.setZero(np, nf);
    out.dz.setZero(np, nf);
  } else {
    out.dx.resize(0, 0);
    out.dy.resize(0, 0);
    out.dz.resize(0, 0);
  }
  for (std::size_t c = 0; c < active.size(); ++c) {
    const auto from = static_cast<Eigen::Index>(c);
    out.value.col(active[c]) = compact.value.col(from);
    if (gradients) {
      out.dx.col(active[c]) = compact.dx.col(from);
      out.dy.col(active[c]) = compact.dy.col(from);
      out.dz.col(active[c]) = compact.dz.col(from);
    }
  }
}

void evaluate_basis_active(const BasisSet& basis, const Vec3* points, std::size_t n_points, bool gradients,
                           BasisValues& out, std::vector<int>& active, double cutoff) {
  const auto np = static_cast<Eigen::Index>(n_points);
  const auto& shells = basis.shells();

  // Shells that contribute in this block, and the column of their first
  // function in the compact matrices.
  std::vector<std::size_t> active_shells;
  std::vector<Eigen::Index> first_column;
  active.clear();
  for (std::size_t s = 0; s < shells.size(); ++s) {
    const Shell& shell = shells[s];
    if (cutoff > 0.0) {
      double r2_min = 1e300;
      for (std::size_t p = 0; p < n_points; ++p) {
        const double dx = points[p][0] - shell.center[0];
        const double dy = points[p][1] - shell.center[1];
        const double dz = points[p][2] - shell.center[2];
        r2_min = std::min(r2_min, dx * dx + dy * dy + dz * dz);
      }
      const Eigen::MatrixXd& T = shell_transformation(shell.l, shell.spherical);
      // |chi| <= sum_c |T_c| |x^a y^b z^c| R <= (largest row sum of |T|) r^l |R|
      const double t_norm = T.cwiseAbs().rowwise().sum().maxCoeff();
      if (t_norm * radial_bound(shell, std::sqrt(r2_min)) < cutoff) continue;
    }
    active_shells.push_back(s);
    first_column.push_back(static_cast<Eigen::Index>(active.size()));
    for (int f = 0; f < shell.n_functions(); ++f) active.push_back(basis.first_function(s) + f);
  }

  const auto nf = static_cast<Eigen::Index>(active.size());
  out.value.resize(np, nf);
  if (gradients) {
    out.dx.resize(np, nf);
    out.dy.resize(np, nf);
    out.dz.resize(np, nf);
  } else {
    out.dx.resize(0, 0);
    out.dy.resize(0, 0);
    out.dz.resize(0, 0);
  }

  const int max_l = basis.max_l();
  // Powers x^0 .. x^(l+1) of the coordinates relative to the shell centre.
  std::vector<double> px(static_cast<std::size_t>(max_l + 2)), py(px.size()), pz(px.size());
  std::vector<double> cart, cart_dx, cart_dy, cart_dz;

  for (std::size_t k = 0; k < active_shells.size(); ++k) {
    const Shell& shell = shells[active_shells[k]];
    const int l = shell.l;

    const auto& powers = cartesian_powers(l);
    const std::size_t nc = powers.size();
    const Eigen::MatrixXd& T = shell_transformation(l, shell.spherical);
    const Eigen::Index first = first_column[k];
    const Eigen::Index n_shell = shell.n_functions();
    cart.resize(nc);
    if (gradients) {
      cart_dx.resize(nc);
      cart_dy.resize(nc);
      cart_dz.resize(nc);
    }

    for (std::size_t p = 0; p < n_points; ++p) {
      const double x = points[p][0] - shell.center[0];
      const double y = points[p][1] - shell.center[1];
      const double z = points[p][2] - shell.center[2];
      const double r2 = x * x + y * y + z * z;

      // Radial part R = sum_k d_k exp(-a_k r^2) and S = sum_k d_k a_k exp(-a_k r^2);
      // dR/dx = -2 x S.
      double R = 0.0, S = 0.0;
      for (std::size_t q = 0; q < shell.n_primitives(); ++q) {
        const double e = shell.norm_coefficients[q] * std::exp(-shell.exponents[q] * r2);
        R += e;
        S += e * shell.exponents[q];
      }

      px[0] = py[0] = pz[0] = 1.0;
      for (int e = 1; e <= l + 1; ++e) {
        const auto i = static_cast<std::size_t>(e);
        px[i] = px[i - 1] * x;
        py[i] = py[i - 1] * y;
        pz[i] = pz[i - 1] * z;
      }

      for (std::size_t c = 0; c < nc; ++c) {
        const auto a = static_cast<std::size_t>(powers[c].x);
        const auto b = static_cast<std::size_t>(powers[c].y);
        const auto d = static_cast<std::size_t>(powers[c].z);
        cart[c] = px[a] * py[b] * pz[d] * R;
        if (gradients) {
          // d/dx [x^a y^b z^d R] = a x^(a-1) y^b z^d R - 2 x^(a+1) y^b z^d S
          cart_dx[c] = (a > 0 ? static_cast<double>(a) * px[a - 1] * R : 0.0) * py[b] * pz[d] -
                       2.0 * px[a + 1] * py[b] * pz[d] * S;
          cart_dy[c] = (b > 0 ? static_cast<double>(b) * py[b - 1] * R : 0.0) * px[a] * pz[d] -
                       2.0 * px[a] * py[b + 1] * pz[d] * S;
          cart_dz[c] = (d > 0 ? static_cast<double>(d) * pz[d - 1] * R : 0.0) * px[a] * py[b] -
                       2.0 * px[a] * py[b] * pz[d + 1] * S;
        }
      }

      const auto row = static_cast<Eigen::Index>(p);
      for (Eigen::Index f = 0; f < n_shell; ++f) {
        double v = 0.0, vx = 0.0, vy = 0.0, vz = 0.0;
        for (std::size_t c = 0; c < nc; ++c) {
          const double t = T(f, static_cast<Eigen::Index>(c));
          if (t == 0.0) continue;
          v += t * cart[c];
          if (gradients) {
            vx += t * cart_dx[c];
            vy += t * cart_dy[c];
            vz += t * cart_dz[c];
          }
        }
        out.value(row, first + f) = v;
        if (gradients) {
          out.dx(row, first + f) = vx;
          out.dy(row, first + f) = vy;
          out.dz(row, first + f) = vz;
        }
      }
    }
  }
}

}  // namespace ellidft
