#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>

#include "core/xyz.hpp"
#include "grid/molecular_grid.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Molecule load(const std::string& stem) {
  return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/" + stem + ".xyz");
}

// Sum of normalised Gaussians (a/pi)^(3/2) exp(-a |r - R_A|^2), one per atom:
// integrates to the number of atoms.
double gaussian_cloud(const Molecule& molecule, const Vec3& r, double a) {
  double value = 0.0;
  for (const Atom& atom : molecule.atoms) {
    const double dx = r[0] - atom.r[0], dy = r[1] - atom.r[1], dz = r[2] - atom.r[2];
    value += std::pow(a / M_PI, 1.5) * std::exp(-a * (dx * dx + dy * dy + dz * dz));
  }
  return value;
}

}  // namespace

TEST_CASE("radial grid integrates Gaussians and exponentials", "[grid][radial]") {
  for (int Z : {1, 6, 8, 17}) {
    const RadialGrid grid = treutler_ahlrichs_grid(75, Z);
    REQUIRE(grid.r.size() == 75);
    for (std::size_t i = 0; i < grid.r.size(); ++i) {
      CHECK(grid.r[i] > 0.0);
      CHECK(grid.weights[i] > 0.0);
      if (i > 0) CHECK(grid.r[i] > grid.r[i - 1]);
    }
    INFO("Z = " << Z);
    // integral r^2 exp(-a r^2) dr = sqrt(pi) / (4 a^(3/2))
    for (double a : {0.3, 1.0, 5.0}) {
      double sum = 0.0;
      for (std::size_t i = 0; i < grid.r.size(); ++i) sum += grid.weights[i] * std::exp(-a * grid.r[i] * grid.r[i]);
      CHECK_THAT(sum, WithinRel(std::sqrt(M_PI) / (4.0 * std::pow(a, 1.5)), 1e-9));
    }
    // integral r^2 exp(-2 r) dr = 2 / 2^3 (density of a hydrogen 1s orbital: has a cusp)
    double sum = 0.0;
    for (std::size_t i = 0; i < grid.r.size(); ++i) sum += grid.weights[i] * std::exp(-2.0 * grid.r[i]);
    CHECK_THAT(sum, WithinRel(0.25, 1e-9));
  }
  CHECK_THROWS_AS(treutler_ahlrichs_grid(0, 1), std::invalid_argument);
  CHECK_THROWS_AS(treutler_ahlrichs_grid(10, 19), std::out_of_range);
  CHECK_THROWS_AS(bragg_radius(0), std::out_of_range);
  CHECK(bragg_radius(6) == 0.70);
}

TEST_CASE("atomic grid: weights and exact integrals", "[grid]") {
  Molecule neon;
  neon.atoms = {{10, {0.4, -0.2, 0.9}}};
  const MolecularGrid grid = build_molecular_grid(neon, {50, 110});
  REQUIRE(grid.size() == 50 * 110);
  double integral = 0.0, dipole_z = 0.0;
  for (std::size_t i = 0; i < grid.size(); ++i) {
    CHECK(grid.weights[i] > 0.0);
    CHECK(grid.atom[i] == 0);
    const double rho = gaussian_cloud(neon, grid.points[i], 0.8);
    integral += grid.weights[i] * rho;
    dipole_z += grid.weights[i] * rho * grid.points[i][2];
  }
  CHECK_THAT(integral, WithinAbs(1.0, 1e-10));
  CHECK_THAT(dipole_z, WithinAbs(0.9, 1e-10));  // centred at the atom
}

TEST_CASE("molecular grid integrates a charge cloud on every atom", "[grid]") {
  for (const char* stem : {"h2o", "hcl", "c2h4"}) {
    const Molecule molecule = load(stem);
    const MolecularGrid grid = build_molecular_grid(molecule, {75, 302});
    INFO(stem);
    REQUIRE(grid.size() == molecule.n_atoms() * 75 * 302);
    double integral = 0.0;
    for (std::size_t i = 0; i < grid.size(); ++i) {
      CHECK(grid.weights[i] >= 0.0);
      integral += grid.weights[i] * gaussian_cloud(molecule, grid.points[i], 1.1);
    }
    // The specification allows 1e-5 for the electron count on this grid
    // (75, 302); the grid is identical to that of the reference code.
    CHECK_THAT(integral, WithinAbs(static_cast<double>(molecule.n_atoms()), 1e-5));
  }
}

TEST_CASE("molecular grid moves rigidly with the molecule", "[grid]") {
  const Molecule molecule = load("h2o");
  Molecule moved = molecule;
  const Vec3 shift = {1.5, -2.5, 0.75};
  for (Atom& atom : moved.atoms)
    for (std::size_t k = 0; k < 3; ++k) atom.r[k] += shift[k];
  // Without ordering by cells, so that point i of one grid is point i of the other.
  const MolecularGrid a = build_molecular_grid(molecule, {20, 50, 0.0});
  const MolecularGrid b = build_molecular_grid(moved, {20, 50, 0.0});
  REQUIRE(a.size() == b.size());
  // Far inside a neighbouring cell the weight is a tiny difference 1 - g with
  // g close to 1, so weights are compared on the scale of the largest weight.
  const double largest = *std::max_element(a.weights.begin(), a.weights.end());
  for (std::size_t i = 0; i < a.size(); ++i) {
    CHECK_THAT(b.weights[i], WithinAbs(a.weights[i], 1e-12 * largest));
    for (std::size_t k = 0; k < 3; ++k) CHECK_THAT(b.points[i][k], WithinAbs(a.points[i][k] + shift[k], 1e-12));
  }
}

TEST_CASE("ordering by cells permutes the points and changes no integral", "[grid]") {
  const Molecule molecule = load("h2o");
  const MolecularGrid plain = build_molecular_grid(molecule, {20, 50, 0.0});
  const MolecularGrid ordered = build_molecular_grid(molecule, {20, 50, 2.5});
  REQUIRE(plain.size() == ordered.size());
  double sum_plain = 0.0, sum_ordered = 0.0, cloud_plain = 0.0, cloud_ordered = 0.0;
  for (std::size_t i = 0; i < plain.size(); ++i) {
    sum_plain += plain.weights[i];
    sum_ordered += ordered.weights[i];
    cloud_plain += plain.weights[i] * gaussian_cloud(molecule, plain.points[i], 0.9);
    cloud_ordered += ordered.weights[i] * gaussian_cloud(molecule, ordered.points[i], 0.9);
  }
  CHECK_THAT(sum_ordered, WithinRel(sum_plain, 1e-12));
  CHECK_THAT(cloud_ordered, WithinRel(cloud_plain, 1e-12));
  // Consecutive points are closer together on average after ordering.
  auto path_length = [](const MolecularGrid& g) {
    double length = 0.0;
    for (std::size_t i = 1; i < g.size(); ++i)
      length += std::hypot(g.points[i][0] - g.points[i - 1][0], g.points[i][1] - g.points[i - 1][1],
                           g.points[i][2] - g.points[i - 1][2]);
    return length;
  };
  CHECK(path_length(ordered) < path_length(plain));
}

TEST_CASE("molecular grid: invalid input", "[grid]") {
  Molecule molecule;
  molecule.atoms = {{1, {0.0, 0.0, 0.0}}, {1, {0.0, 0.0, 0.0}}};
  CHECK_THROWS_AS(build_molecular_grid(molecule, {10, 14}), std::invalid_argument);
  molecule.atoms[1].r[2] = 1.4;
  CHECK_THROWS_AS(build_molecular_grid(molecule, {10, 15}), std::invalid_argument);  // no such Lebedev grid
  CHECK(build_molecular_grid(molecule, {10, 14}).size() == 280);
}
