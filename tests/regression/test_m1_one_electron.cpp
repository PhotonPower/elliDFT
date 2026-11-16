// Milestone M1: overlap, kinetic, nuclear attraction and dipole matrices
// against the reference code. Acceptance criterion: largest deviation at
// most 1e-10.
//
// Reference: tests/reference/data/m1/ints1e_*.json and highl_*.json
// (PySCF, tests/reference/gen_m1.py).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "integrals/one_electron.hpp"
#include "reference_data.hpp"

using namespace ellidft;
using nlohmann::json;

namespace {

constexpr double acceptance_threshold = 1e-10;

struct Deviations {
  double overlap = 0.0, kinetic = 0.0, nuclear = 0.0, dipole = 0.0;
  double largest() const { return std::max({overlap, kinetic, nuclear, dipole}); }
};

double max_deviation(const Eigen::MatrixXd& computed, const json& packed) {
  const Eigen::MatrixXd reference = testing::unpack_symmetric(packed, computed.rows());
  return (computed - reference).cwiseAbs().maxCoeff();
}

Deviations compare(const json& record) {
  const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), "reference");
  const BasisSet basis = testing::basis_from_reference(record, molecule);
  REQUIRE(basis.n_functions() == record.at("n_functions").get<int>());
  const json& matrices = record.at("matrices");

  Deviations d;
  d.overlap = max_deviation(overlap_matrix(basis), matrices.at("overlap"));
  d.kinetic = max_deviation(kinetic_matrix(basis), matrices.at("kinetic"));
  d.nuclear = max_deviation(nuclear_attraction_matrix(basis, molecule), matrices.at("nuclear"));
  const auto dipole = dipole_matrices(basis);
  const char* const keys[3] = {"dipole_x", "dipole_y", "dipole_z"};
  for (std::size_t k = 0; k < 3; ++k) d.dipole = std::max(d.dipole, max_deviation(dipole[k], matrices.at(keys[k])));
  return d;
}

}  // namespace

TEST_CASE("M1: one-electron matrices of the test molecules equal the reference", "[m1][regression]") {
  const auto files = testing::reference_files("m1", "ints1e_");
  REQUIRE(files.size() == 57);
  Deviations worst;
  for (const std::string& name : files) {
    const Deviations d = compare(testing::load_reference("m1/" + name));
    INFO(name << ": overlap " << d.overlap << ", kinetic " << d.kinetic << ", nuclear " << d.nuclear
              << ", dipole " << d.dipole);
    CHECK(d.overlap <= acceptance_threshold);
    CHECK(d.kinetic <= acceptance_threshold);
    CHECK(d.nuclear <= acceptance_threshold);
    CHECK(d.dipole <= acceptance_threshold);
    worst.overlap = std::max(worst.overlap, d.overlap);
    worst.kinetic = std::max(worst.kinetic, d.kinetic);
    worst.nuclear = std::max(worst.nuclear, d.nuclear);
    worst.dipole = std::max(worst.dipole, d.dipole);
  }
  // Shown with the option -s; quoted in the milestone report.
  SUCCEED("largest deviations over 57 cases: overlap " << worst.overlap << ", kinetic " << worst.kinetic
                                                      << ", nuclear " << worst.nuclear << ", dipole "
                                                      << worst.dipole);
}

// Artificial basis with angular momenta 0..6 on three centres without
// symmetry: every pair of angular momenta occurs between different centres.
TEST_CASE("M1: one-electron matrices up to l = 6 equal the reference", "[m1][regression]") {
  for (const std::string name : {"highl_sph.json", "highl_cart.json"}) {
    DYNAMIC_SECTION(name) {
      const Deviations d = compare(testing::load_reference("m1/" + name));
      INFO("overlap " << d.overlap << ", kinetic " << d.kinetic << ", nuclear " << d.nuclear << ", dipole "
                      << d.dipole);
      CHECK(d.overlap <= acceptance_threshold);
      CHECK(d.kinetic <= acceptance_threshold);
      CHECK(d.nuclear <= acceptance_threshold);
      CHECK(d.dipole <= acceptance_threshold);
      SUCCEED(name << ": overlap " << d.overlap << ", kinetic " << d.kinetic << ", nuclear " << d.nuclear
                   << ", dipole " << d.dipole);
    }
  }
}
