#pragma once

#include <array>
#include <vector>

namespace ellidft {

// A quadrature on the unit sphere: integral f dOmega / (4 pi) ~ sum_i w_i f(p_i).
// The weights sum to 1.
struct AngularGrid {
  std::vector<std::array<double, 3>> points;  // unit vectors
  std::vector<double> weights;
  int degree = 0;  // polynomials up to this degree are integrated exactly
};

// Available Lebedev grids as numbers of points: 6, 14, 26, ..., 5810.
const std::vector<int>& lebedev_sizes();

// Nominal algebraic degree of the Lebedev grid with `n_points` points.
// Throws std::invalid_argument if there is no such grid.
int lebedev_degree(int n_points);

// Lebedev grid with `n_points` points, from the tables of the Lebedev-Laikov
// routine (third_party/lebedev). Throws std::invalid_argument if there is no
// such grid.
AngularGrid lebedev_grid(int n_points);

}  // namespace ellidft
