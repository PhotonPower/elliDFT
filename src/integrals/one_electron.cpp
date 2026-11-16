#include "integrals/one_electron.hpp"

#include <cmath>

#include "basis/solid_harmonics.hpp"
#include "integrals/hermite.hpp"
#include "integrals/shell_pair.hpp"

namespace ellidft {

namespace {

// Everything a kernel needs for one pair of primitives. The Hermite
// coefficients are those of the three cartesian directions; only t = 0
// enters integrals over all space (Hernández Vera 2025, Eq. 15).
struct PrimitiveContext {
  const HermiteE& Ex;
  const HermiteE& Ey;
  const HermiteE& Ez;
  double p;       // sum of the exponents
  double b;       // exponent of the second primitive
  const Vec3& B;  // centre of the second shell
};

// Integrals over one shell pair for N operators at once. The kernel fills
// value[0..N) for one pair of cartesian Gaussians (x-A)^pa exp(-a r_A^2),
// (x-B)^pb exp(-b r_B^2); contraction and the transformation to the basis
// functions are done here. `extra_j` is how far the kernel raises the power
// of the second function.
template <std::size_t N, typename Kernel>
std::array<Eigen::MatrixXd, N> shell_block(const Shell& a, const Shell& b, int extra_j, Kernel&& kernel) {
  const auto& powers_a = cartesian_powers(a.l);
  const auto& powers_b = cartesian_powers(b.l);
  const auto nca = static_cast<Eigen::Index>(powers_a.size());
  const auto ncb = static_cast<Eigen::Index>(powers_b.size());

  std::array<Eigen::MatrixXd, N> raw;
  for (auto& m : raw) m = Eigen::MatrixXd::Zero(nca, ncb);

  for (std::size_t ka = 0; ka < a.n_primitives(); ++ka) {
    for (std::size_t kb = 0; kb < b.n_primitives(); ++kb) {
      const double ea = a.exponents[ka];
      const double eb = b.exponents[kb];
      const double coefficient = a.norm_coefficients[ka] * b.norm_coefficients[kb];
      const HermiteE Ex(a.l, b.l + extra_j, ea, eb, a.center[0], b.center[0]);
      const HermiteE Ey(a.l, b.l + extra_j, ea, eb, a.center[1], b.center[1]);
      const HermiteE Ez(a.l, b.l + extra_j, ea, eb, a.center[2], b.center[2]);
      const PrimitiveContext context{Ex, Ey, Ez, ea + eb, eb, b.center};

      std::array<double, N> value{};
      for (Eigen::Index ca = 0; ca < nca; ++ca)
        for (Eigen::Index cb = 0; cb < ncb; ++cb) {
          kernel(context, powers_a[static_cast<std::size_t>(ca)], powers_b[static_cast<std::size_t>(cb)], value);
          for (std::size_t k = 0; k < N; ++k) raw[k](ca, cb) += coefficient * value[k];
        }
    }
  }

  const Eigen::MatrixXd& Ta = shell_transformation(a.l, a.spherical);
  const Eigen::MatrixXd& Tb = shell_transformation(b.l, b.spherical);
  std::array<Eigen::MatrixXd, N> out;
  for (std::size_t k = 0; k < N; ++k) out[k] = Ta * raw[k] * Tb.transpose();
  return out;
}

template <std::size_t N, typename Kernel>
std::array<Eigen::MatrixXd, N> assemble(const BasisSet& basis, int extra_j, Kernel&& kernel) {
  const Eigen::Index n = basis.n_functions();
  std::array<Eigen::MatrixXd, N> matrices;
  for (auto& m : matrices) m = Eigen::MatrixXd::Zero(n, n);
  const auto& shells = basis.shells();
  for (std::size_t sa = 0; sa < shells.size(); ++sa)
    for (std::size_t sb = 0; sb < shells.size(); ++sb) {
      const auto blocks = shell_block<N>(shells[sa], shells[sb], extra_j, kernel);
      for (std::size_t k = 0; k < N; ++k)
        matrices[k].block(basis.first_function(sa), basis.first_function(sb), blocks[k].rows(),
                          blocks[k].cols()) = blocks[k];
    }
  return matrices;
}

// (pi/p)^(3/2): product of the three one-dimensional factors sqrt(pi/p)
// (Hernández Vera 2025, Eqs. 35, 36).
double volume_factor(double p) { return std::pow(M_PI / p, 1.5); }

// One-dimensional kinetic energy integral without the factor sqrt(pi/p)
// (Hernández Vera 2025, Eq. 42):
//   -2 b^2 E_0^{i,j+2} + b (2j+1) E_0^{ij} - j (j-1)/2 E_0^{i,j-2}
double kinetic_1d(const HermiteE& E, int i, int j, double b) {
  double value = -2.0 * b * b * E(i, j + 2, 0) + b * (2 * j + 1) * E(i, j, 0);
  if (j >= 2) value -= 0.5 * j * (j - 1) * E(i, j - 2, 0);
  return value;
}

}  // namespace

Eigen::MatrixXd overlap_matrix(const BasisSet& basis) {
  // Hernández Vera 2025, Eq. 36.
  auto kernel = [](const PrimitiveContext& c, const CartesianPowers& pa, const CartesianPowers& pb,
                   std::array<double, 1>& value) {
    value[0] = volume_factor(c.p) * c.Ex(pa.x, pb.x, 0) * c.Ey(pa.y, pb.y, 0) * c.Ez(pa.z, pb.z, 0);
  };
  return assemble<1>(basis, 0, kernel)[0];
}

Eigen::MatrixXd kinetic_matrix(const BasisSet& basis) {
  // Hernández Vera 2025, Eqs. 38 and 43: T = T_x S_y S_z + S_x T_y S_z + S_x S_y T_z.
  auto kernel = [](const PrimitiveContext& c, const CartesianPowers& pa, const CartesianPowers& pb,
                   std::array<double, 1>& value) {
    const double sx = c.Ex(pa.x, pb.x, 0), sy = c.Ey(pa.y, pb.y, 0), sz = c.Ez(pa.z, pb.z, 0);
    const double tx = kinetic_1d(c.Ex, pa.x, pb.x, c.b);
    const double ty = kinetic_1d(c.Ey, pa.y, pb.y, c.b);
    const double tz = kinetic_1d(c.Ez, pa.z, pb.z, c.b);
    value[0] = volume_factor(c.p) * (tx * sy * sz + sx * ty * sz + sx * sy * tz);
  };
  return assemble<1>(basis, 2, kernel)[0];
}

std::array<Eigen::MatrixXd, 3> dipole_matrices(const BasisSet& basis, const Vec3& origin) {
  // With x - C = (x - B) + (B - C), the factor (x - B) raises the power of
  // the second function by one:
  //   <i| x - C |j> = S_{i,j+1} + (B - C) S_{ij}
  // in the direction of the operator, times the overlaps of the other two.
  auto kernel = [&origin](const PrimitiveContext& c, const CartesianPowers& pa, const CartesianPowers& pb,
                          std::array<double, 3>& value) {
    const double sx = c.Ex(pa.x, pb.x, 0), sy = c.Ey(pa.y, pb.y, 0), sz = c.Ez(pa.z, pb.z, 0);
    const double dx = c.Ex(pa.x, pb.x + 1, 0) + (c.B[0] - origin[0]) * sx;
    const double dy = c.Ey(pa.y, pb.y + 1, 0) + (c.B[1] - origin[1]) * sy;
    const double dz = c.Ez(pa.z, pb.z + 1, 0) + (c.B[2] - origin[2]) * sz;
    const double f = volume_factor(c.p);
    value[0] = f * dx * sy * sz;
    value[1] = f * sx * dy * sz;
    value[2] = f * sx * sy * dz;
  };
  return assemble<3>(basis, 1, kernel);
}

Eigen::MatrixXd nuclear_attraction_matrix(const BasisSet& basis, const Molecule& molecule) {
  const Eigen::Index n = basis.n_functions();
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(n, n);
  const auto& shells = basis.shells();
  HermiteCoulomb R(2 * basis.max_l());

  for (std::size_t sa = 0; sa < shells.size(); ++sa)
    for (std::size_t sb = 0; sb < shells.size(); ++sb) {
      const ShellPair pair = make_shell_pair(shells[sa], shells[sb]);
      const std::vector<HermiteIndex> hermite = hermite_indices(pair.L);
      const auto nf = static_cast<std::size_t>(pair.nf);
      std::vector<double> block(nf, 0.0);

      for (const PrimitivePair& prim : pair.primitives) {
        for (const Atom& atom : molecule.atoms) {
          // Hernández Vera 2025, Eq. 49:
          //   <a| 1/r_C |b> = (2 pi / p) sum_tuv E_tuv R^0_tuv(p, P - C)
          R.compute(pair.L, prim.p, prim.P[0] - atom.r[0], prim.P[1] - atom.r[1], prim.P[2] - atom.r[2]);
          const double factor = -atom.Z * 2.0 * M_PI / prim.p;
          for (std::size_t h = 0; h < hermite.size(); ++h) {
            const double r = factor * R.value(hermite[h].t, hermite[h].u, hermite[h].v);
            const double* E = &prim.E[h * nf];
            for (std::size_t f = 0; f < nf; ++f) block[f] += r * E[f];
          }
        }
      }
      const int ia = basis.first_function(sa), ib = basis.first_function(sb);
      for (int fa = 0; fa < pair.nfa; ++fa)
        for (int fb = 0; fb < pair.nfb; ++fb)
          V(ia + fa, ib + fb) = block[static_cast<std::size_t>(fa * pair.nfb + fb)];
    }
  return V;
}

Eigen::MatrixXd overlap_matrix_from_shell_pairs(const BasisSet& basis) {
  const Eigen::Index n = basis.n_functions();
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  const auto& shells = basis.shells();
  for (std::size_t sa = 0; sa < shells.size(); ++sa)
    for (std::size_t sb = 0; sb < shells.size(); ++sb) {
      const ShellPair pair = make_shell_pair(shells[sa], shells[sb]);
      const int ia = basis.first_function(sa), ib = basis.first_function(sb);
      // Only the Hermite Gaussian with t = u = v = 0 has a non-zero integral,
      // (pi/p)^(3/2); it is the first entry of hermite_indices().
      for (const PrimitivePair& prim : pair.primitives)
        for (int fa = 0; fa < pair.nfa; ++fa)
          for (int fb = 0; fb < pair.nfb; ++fb)
            S(ia + fa, ib + fb) += volume_factor(prim.p) * prim.E[static_cast<std::size_t>(fa * pair.nfb + fb)];
    }
  return S;
}

}  // namespace ellidft
