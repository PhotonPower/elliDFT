// Milestone M2: basis functions and their gradients at points against the
// reference code.
//
// Reference: tests/reference/data/m2/basis_values.json (PySCF eval_gto).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "grid/basis_on_grid.hpp"
#include "reference_data.hpp"

using namespace ellidft;
using nlohmann::json;

namespace {

double max_deviation(const Eigen::MatrixXd& computed, const json& reference) {
  double worst = 0.0;
  for (Eigen::Index p = 0; p < computed.rows(); ++p)
    for (Eigen::Index f = 0; f < computed.cols(); ++f)
      worst = std::max(worst, std::abs(computed(p, f) - reference[static_cast<std::size_t>(p)]
                                                                 [static_cast<std::size_t>(f)].get<double>()));
  return worst;
}

}  // namespace

TEST_CASE("M2: basis functions and gradients at points equal the reference", "[m2][regression][grid]") {
  const json reference = testing::load_reference("m2/basis_values.json");
  std::vector<Vec3> points;
  for (const json& p : reference.at("points_bohr")) points.push_back({p[0].get<double>(), p[1].get<double>(), p[2].get<double>()});
  REQUIRE(reference.at("cases").size() == 4);

  for (const json& record : reference.at("cases")) {
    const std::string label = record.at("basis").at("name").get<std::string>() +
                              (record.at("basis").at("spherical").get<bool>() ? " spherical" : " cartesian");
    DYNAMIC_SECTION(label) {
      const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), label);
      const BasisSet basis = testing::basis_from_reference(record, molecule);
      REQUIRE(basis.n_functions() == record.at("n_functions").get<int>());

      BasisValues values;
      evaluate_basis(basis, points.data(), points.size(), true, values, 0.0);
      const double dv = max_deviation(values.value, record.at("value"));
      const double dx = max_deviation(values.dx, record.at("dx"));
      const double dy = max_deviation(values.dy, record.at("dy"));
      const double dz = max_deviation(values.dz, record.at("dz"));
      INFO("value " << dv << ", gradient " << dx << " " << dy << " " << dz);
      CHECK(dv < 1e-11);
      CHECK(dx < 1e-10);
      CHECK(dy < 1e-10);
      CHECK(dz < 1e-10);
      SUCCEED(label << ": value " << dv << ", gradient " << std::max({dx, dy, dz}));

      // With the default cutoff only negligible values are set to zero.
      BasisValues screened;
      evaluate_basis(basis, points.data(), points.size(), true, screened);
      CHECK((screened.value - values.value).cwiseAbs().maxCoeff() < 1e-15);
      CHECK((screened.dx - values.dx).cwiseAbs().maxCoeff() < 1e-13);
    }
  }
}
