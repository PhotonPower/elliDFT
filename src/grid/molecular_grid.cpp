#include "grid/molecular_grid.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include "core/elements.hpp"
#include "grid/lebedev.hpp"

// The formulas of this file follow the grid construction of the reference
// code, read in the PySCF 2.15.0 sources dft/radi.py and dft/gen_grid.py
// (see docs/sources/README.md). PySCF attributes them to Treutler and
// Ahlrichs, J. Chem. Phys. 102, 346 (1995), and Becke, J. Chem. Phys. 88,
// 2547 (1988); those papers were not consulted.

namespace ellidft {

namespace {

void check_element(int Z) {
  if (Z < 1 || Z > max_atomic_number)
    throw std::out_of_range("grid: unsupported atomic number " + std::to_string(Z));
}

// Radius parameter xi of the Treutler-Ahlrichs radial grid, H to Ar
// (PySCF dft/radi.py, _treutler_ahlrichs_xi).
constexpr double treutler_xi[max_atomic_number] = {0.8, 0.9, 1.8, 1.4, 1.3, 1.1, 0.9, 0.9, 0.9,
                                                   0.9, 1.4, 1.3, 1.3, 1.2, 1.1, 1.0, 1.0, 1.0};

// Bragg-Slater radii in Angstrom, H to Ar (PySCF data/radii.py, BRAGG).
constexpr double bragg[max_atomic_number] = {0.35, 1.40, 1.45, 1.05, 0.85, 0.70, 0.65, 0.60, 0.50,
                                             1.50, 1.80, 1.50, 1.25, 1.10, 1.00, 1.00, 1.00, 1.80};

// Becke's cell function: the polynomial g -> g (3 - g^2) / 2 applied three
// times (PySCF dft/gen_grid.py, original_becke).
double becke_step(double g) {
  for (int k = 0; k < 3; ++k) g = 0.5 * g * (3.0 - g * g);
  return g;
}

}  // namespace

double bragg_radius(int Z) {
  check_element(Z);
  return bragg[Z - 1];
}

RadialGrid treutler_ahlrichs_grid(int n_points, int Z) {
  check_element(Z);
  if (n_points < 1) throw std::invalid_argument("treutler_ahlrichs_grid: at least one point is needed");
  RadialGrid grid;
  grid.r.resize(static_cast<std::size_t>(n_points));
  grid.weights.resize(static_cast<std::size_t>(n_points));
  // PySCF dft/radi.py, treutler_ahlrichs: x_i = cos(i pi / (n+1)), i = 1..n,
  //   r = -xi/ln2 (1+x)^0.6 ln((1-x)/2),
  //   dr = pi/(n+1) sin(i pi/(n+1)) xi/ln2 (1+x)^0.6 [-0.6/(1+x) ln((1-x)/2) + 1/(1-x)].
  // The points are stored with increasing r.
  const double step = M_PI / (n_points + 1);
  const double scale = treutler_xi[Z - 1] / std::log(2.0);
  for (int i = 0; i < n_points; ++i) {
    const double x = std::cos((i + 1) * step);
    const double power = std::pow(1.0 + x, 0.6);
    const double log_term = std::log(0.5 * (1.0 - x));
    const double r = -scale * power * log_term;
    const double dr = step * std::sin((i + 1) * step) * scale * power * (-0.6 / (1.0 + x) * log_term + 1.0 / (1.0 - x));
    const auto k = static_cast<std::size_t>(n_points - 1 - i);
    grid.r[k] = r;
    grid.weights[k] = r * r * dr;
  }
  return grid;
}

MolecularGrid build_molecular_grid(const Molecule& molecule, const GridSpec& spec) {
  const std::size_t n_atoms = molecule.n_atoms();
  const AngularGrid angular = lebedev_grid(spec.n_angular);

  // Size adjustment (PySCF dft/radi.py, treutler_atomic_radii_adjust):
  //   a_ij = (rho_j/rho_i - rho_i/rho_j) / 4 with rho = sqrt(Bragg radius),
  //   limited to [-1/2, 1/2];  g -> g + a_ij (1 - g^2).
  std::vector<double> adjust(n_atoms * n_atoms, 0.0), inverse_distance(n_atoms * n_atoms, 0.0);
  for (std::size_t i = 0; i < n_atoms; ++i)
    for (std::size_t j = 0; j < n_atoms; ++j) {
      if (i == j) continue;
      const double ri = std::sqrt(bragg_radius(molecule.atoms[i].Z));
      const double rj = std::sqrt(bragg_radius(molecule.atoms[j].Z));
      adjust[i * n_atoms + j] = std::clamp(0.25 * (rj / ri - ri / rj), -0.5, 0.5);
      const double dx = molecule.atoms[i].r[0] - molecule.atoms[j].r[0];
      const double dy = molecule.atoms[i].r[1] - molecule.atoms[j].r[1];
      const double dz = molecule.atoms[i].r[2] - molecule.atoms[j].r[2];
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (distance == 0.0) throw std::invalid_argument("build_molecular_grid: two atoms at the same position");
      inverse_distance[i * n_atoms + j] = 1.0 / distance;
    }

  MolecularGrid grid;
  const std::size_t per_atom = static_cast<std::size_t>(spec.n_radial) * angular.points.size();
  grid.points.reserve(n_atoms * per_atom);
  grid.weights.reserve(n_atoms * per_atom);
  grid.atom.reserve(n_atoms * per_atom);

  std::vector<double> distance(n_atoms), cell(n_atoms);
  for (std::size_t a = 0; a < n_atoms; ++a) {
    const RadialGrid radial = treutler_ahlrichs_grid(spec.n_radial, molecule.atoms[a].Z);
    for (std::size_t ir = 0; ir < radial.r.size(); ++ir) {
      for (std::size_t ia = 0; ia < angular.points.size(); ++ia) {
        Vec3 point;
        for (std::size_t k = 0; k < 3; ++k)
          point[k] = molecule.atoms[a].r[k] + radial.r[ir] * angular.points[ia][k];

        // Becke partitioning (PySCF dft/gen_grid.py, get_partition): for
        // every pair i > j, g = (|r - R_i| - |r - R_j|) / R_ij, adjusted and
        // smoothed; P_i *= (1 - g)/2, P_j *= (1 + g)/2. The weight of atom a
        // is P_a / sum_k P_k.
        for (std::size_t k = 0; k < n_atoms; ++k) {
          const double dx = point[0] - molecule.atoms[k].r[0];
          const double dy = point[1] - molecule.atoms[k].r[1];
          const double dz = point[2] - molecule.atoms[k].r[2];
          distance[k] = std::sqrt(dx * dx + dy * dy + dz * dz);
          cell[k] = 1.0;
        }
        for (std::size_t i = 0; i < n_atoms; ++i)
          for (std::size_t j = 0; j < i; ++j) {
            double g = (distance[i] - distance[j]) * inverse_distance[i * n_atoms + j];
            g += adjust[i * n_atoms + j] * (1.0 - g * g);
            g = becke_step(g);
            cell[i] *= 0.5 * (1.0 - g);
            cell[j] *= 0.5 * (1.0 + g);
          }
        double total = 0.0;
        for (std::size_t k = 0; k < n_atoms; ++k) total += cell[k];

        // Solid angle: the Lebedev weights sum to 1, hence the factor 4 pi.
        const double volume = 4.0 * M_PI * radial.weights[ir] * angular.weights[ia];
        grid.points.push_back(point);
        grid.weights.push_back(volume * cell[a] / total);
        grid.atom.push_back(static_cast<int>(a));
      }
    }
  }

  if (spec.cell_size > 0.0 && grid.size() > 0) {
    // Sort by cell index (x, then y, then z); stable, so the order inside a
    // cell stays atom by atom.
    Vec3 lower = grid.points[0];
    for (const Vec3& p : grid.points)
      for (std::size_t k = 0; k < 3; ++k) lower[k] = std::min(lower[k], p[k]);
    auto cell_index = [&](const Vec3& p) {
      std::array<long, 3> c{};
      for (std::size_t k = 0; k < 3; ++k) c[k] = static_cast<long>((p[k] - lower[k]) / spec.cell_size);
      return c;
    };
    std::vector<std::size_t> order(grid.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return cell_index(grid.points[a]) < cell_index(grid.points[b]); });
    MolecularGrid sorted;
    sorted.points.reserve(grid.size());
    sorted.weights.reserve(grid.size());
    sorted.atom.reserve(grid.size());
    for (std::size_t i : order) {
      sorted.points.push_back(grid.points[i]);
      sorted.weights.push_back(grid.weights[i]);
      sorted.atom.push_back(grid.atom[i]);
    }
    return sorted;
  }
  return grid;
}

}  // namespace ellidft
