#pragma once

#include <vector>

#include "core/molecule.hpp"

namespace ellidft {

// Size of the atomic grids: radial points and Lebedev points per radial shell.
// No pruning: every radial shell carries the full angular grid.
struct GridSpec {
  int n_radial = 75;
  int n_angular = 302;
  // The points are ordered by cubic cells of this edge length (Bohr), so that
  // consecutive points are close in space and blocks of points see few basis
  // functions. 0 keeps the order atom by atom, shell by shell. The order does
  // not change any integral beyond rounding.
  double cell_size = 2.5;
};

// Radial quadrature of Treutler and Ahlrichs (mapping M4 with Chebyshev
// points of the second kind), as implemented in the reference code:
//
//   integral_0^inf f(r) r^2 dr  ~  sum_i weights[i] f(r[i]).
//
// The scale parameter xi depends on the element (H to Ar).
struct RadialGrid {
  std::vector<double> r;
  std::vector<double> weights;  // include the factor r^2
};
RadialGrid treutler_ahlrichs_grid(int n_points, int Z);

// Bragg-Slater radius of an element in Angstrom (H to Ar), used for the size
// adjustment of the cell boundaries. Only ratios of these radii enter.
double bragg_radius(int Z);

// Integration grid of a molecule: integral f(r) d^3r ~ sum_i weights[i] f(points[i]).
struct MolecularGrid {
  std::vector<Vec3> points;      // Bohr
  std::vector<double> weights;
  std::vector<int> atom;         // atom whose atomic grid the point belongs to
  std::size_t size() const { return points.size(); }
};

// Atomic grids (radial times Lebedev) combined with Becke's partitioning of
// space into fuzzy atomic cells, with the size adjustment of Treutler and
// Ahlrichs. Same construction as the default of the reference code with
// pruning switched off.
MolecularGrid build_molecular_grid(const Molecule& molecule, const GridSpec& spec = {});

}  // namespace ellidft
