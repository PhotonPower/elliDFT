#include "grid/lebedev.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace ellidft {

namespace {

struct LebedevOrbit {
  int code;
  double a, b, weight;
};

struct LebedevTable {
  int n_points;
  int n_orbits;
  const LebedevOrbit* orbits;
};

#include "grid/lebedev_tables.inc"

// Nominal algebraic degrees of the grids, in the order of lebedev_tables.
// Values as tabulated in PySCF 2.15.0, dft/LebedevGrid.py (LEBEDEV_ORDER),
// consulted 2026-11-14; the unit tests verify them by integrating all
// monomials up to the degree.
constexpr int degrees[] = {3,  5,  7,  9,  11, 13, 15, 17, 19, 21,  23,  25,  27,  29,  31,  35,
                           41, 47, 53, 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131};

constexpr std::size_t n_tables = sizeof(lebedev_tables) / sizeof(lebedev_tables[0]);
static_assert(n_tables == sizeof(degrees) / sizeof(degrees[0]));

std::size_t table_index(int n_points) {
  for (std::size_t i = 0; i < n_tables; ++i)
    if (lebedev_tables[i].n_points == n_points) return i;
  throw std::invalid_argument("no Lebedev grid with " + std::to_string(n_points) + " points");
}

// All points equivalent to (x, y, z) under the octahedral group with
// inversion: the six permutations of the coordinates with all sign
// combinations, without duplicates. This reproduces the orbits listed in the
// header of GEN_OH in Lebedev-Laikov.F:
//   code 1: (0,0,1)                      6 points
//   code 2: (0,a,a), a = 1/sqrt(2)      12 points
//   code 3: (a,a,a), a = 1/sqrt(3)       8 points
//   code 4: (a,a,b), b = sqrt(1 - 2a^2) 24 points
//   code 5: (a,b,0), b = sqrt(1 - a^2)  24 points
//   code 6: (a,b,c), c = sqrt(1-a^2-b^2) 48 points
std::vector<std::array<double, 3>> octahedral_orbit(double x, double y, double z) {
  std::vector<std::array<double, 3>> points;
  std::array<double, 3> base = {x, y, z};
  std::sort(base.begin(), base.end());
  do {
    for (int signs = 0; signs < 8; ++signs) {
      std::array<double, 3> p = base;
      for (std::size_t k = 0; k < 3; ++k)
        if (signs & (1 << k)) p[k] = -p[k];
      for (double& c : p) c += 0.0;  // -0.0 -> 0.0
      points.push_back(p);
    }
  } while (std::next_permutation(base.begin(), base.end()));
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  return points;
}

}  // namespace

const std::vector<int>& lebedev_sizes() {
  static const std::vector<int> sizes = [] {
    std::vector<int> s;
    for (std::size_t i = 0; i < n_tables; ++i) s.push_back(lebedev_tables[i].n_points);
    return s;
  }();
  return sizes;
}

int lebedev_degree(int n_points) { return degrees[table_index(n_points)]; }

AngularGrid lebedev_grid(int n_points) {
  const std::size_t index = table_index(n_points);
  const LebedevTable& table = lebedev_tables[index];
  AngularGrid grid;
  grid.degree = degrees[index];
  grid.points.reserve(static_cast<std::size_t>(n_points));
  grid.weights.reserve(static_cast<std::size_t>(n_points));

  for (int i = 0; i < table.n_orbits; ++i) {
    const LebedevOrbit& orbit = table.orbits[i];
    double x = 0.0, y = 0.0, z = 0.0;
    std::size_t expected = 0;
    switch (orbit.code) {
      case 1: x = 0.0; y = 0.0; z = 1.0; expected = 6; break;
      case 2: x = 0.0; y = z = std::sqrt(0.5); expected = 12; break;
      case 3: x = y = z = std::sqrt(1.0 / 3.0); expected = 8; break;
      case 4: x = y = orbit.a; z = std::sqrt(1.0 - 2.0 * orbit.a * orbit.a); expected = 24; break;
      case 5: x = orbit.a; y = std::sqrt(1.0 - orbit.a * orbit.a); z = 0.0; expected = 24; break;
      case 6:
        x = orbit.a; y = orbit.b; z = std::sqrt(1.0 - orbit.a * orbit.a - orbit.b * orbit.b);
        expected = 48;
        break;
      default: throw std::logic_error("lebedev_grid: invalid orbit code in the tables");
    }
    const auto points = octahedral_orbit(x, y, z);
    if (points.size() != expected) throw std::logic_error("lebedev_grid: orbit has the wrong number of points");
    for (const auto& p : points) {
      grid.points.push_back(p);
      grid.weights.push_back(orbit.weight);
    }
  }
  if (grid.points.size() != static_cast<std::size_t>(n_points))
    throw std::logic_error("lebedev_grid: wrong number of points");
  return grid;
}

}  // namespace ellidft
