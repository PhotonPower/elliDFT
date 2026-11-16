// Milestone M2: molecular grids point by point against the reference code.
//
// Reference: tests/reference/data/m2/grid_*.json (PySCF, tests/reference/gen_m2.py).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "grid/molecular_grid.hpp"
#include "reference_data.hpp"

using namespace ellidft;
using nlohmann::json;

TEST_CASE("M2: molecular grids equal the reference point by point", "[m2][regression][grid]") {
  const auto files = testing::reference_files("m2", "grid_");
  REQUIRE(files.size() == 3);
  for (const std::string& name : files) {
    DYNAMIC_SECTION(name) {
      const json record = testing::load_reference("m2/" + name);
      const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), name);
      const GridSpec spec{record.at("n_radial").get<int>(), record.at("n_angular").get<int>()};
      const MolecularGrid grid = build_molecular_grid(molecule, spec);

      const auto points = record.at("points").get<std::vector<std::vector<double>>>();
      const auto weights = record.at("weights").get<std::vector<double>>();
      const auto atoms = record.at("atom").get<std::vector<int>>();
      REQUIRE(grid.size() == points.size());

      // The two programs order the points differently. Every reference point
      // is matched with the nearest point of the same atom.
      std::vector<bool> used(grid.size(), false);
      double worst_position = 0.0, worst_weight = 0.0, largest_weight = 0.0;
      for (std::size_t r = 0; r < points.size(); ++r) {
        std::size_t best = grid.size();
        double best_distance = 1e300;
        for (std::size_t i = 0; i < grid.size(); ++i) {
          if (used[i] || grid.atom[i] != atoms[r]) continue;
          const double d = std::hypot(grid.points[i][0] - points[r][0], grid.points[i][1] - points[r][1],
                                      grid.points[i][2] - points[r][2]);
          if (d < best_distance) {
            best_distance = d;
            best = i;
          }
        }
        REQUIRE(best < grid.size());
        used[best] = true;
        worst_position = std::max(worst_position, best_distance);
        worst_weight = std::max(worst_weight, std::abs(grid.weights[best] - weights[r]));
        largest_weight = std::max(largest_weight, std::abs(weights[r]));
      }
      INFO("largest position deviation " << worst_position << ", largest weight deviation " << worst_weight
                                         << " (largest weight " << largest_weight << ")");
      CHECK(worst_position < 1e-11);
      CHECK(worst_weight < 1e-12 * largest_weight);
      SUCCEED(name << ": position " << worst_position << ", weight " << worst_weight << " of " << largest_weight);
    }
  }
}
