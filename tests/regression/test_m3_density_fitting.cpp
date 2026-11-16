// Milestone M3: density fitting against the reference code.
// Acceptance criterion: energies with density fitting within 1e-7 Hartree of
// the reference with the same auxiliary basis.
//
// Reference: tests/reference/data/m3/*.json (PySCF, tests/reference/gen_m3.py).
// Every test run covers the atoms and H2; the other cases are tagged [slow].

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "df/df_integrals.hpp"
#include "df/df_jk.hpp"
#include "reference_data.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using nlohmann::json;

TEST_CASE("M3: two- and three-centre integrals equal the reference", "[m3][regression][df]") {
  const auto files = testing::reference_files("m3", "df_ints_");
  REQUIRE(files.size() == 2);
  for (const std::string& name : files) {
    DYNAMIC_SECTION(name) {
      const json record = testing::load_reference("m3/" + name);
      const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), name);
      const BasisSet basis = testing::basis_from_reference(record, molecule);
      const BasisSet aux =
          BasisSet::build(molecule, BasisLibrary::load(record.at("auxiliary_basis").get<std::string>()), true);
      REQUIRE(basis.n_functions() == record.at("n_functions").get<int>());
      REQUIRE(aux.n_functions() == record.at("n_auxiliary").get<int>());

      const Eigen::MatrixXd V = two_center_coulomb(aux);
      const double two = (V - testing::unpack_symmetric(record.at("two_center"), V.rows())).cwiseAbs().maxCoeff();

      const Eigen::MatrixXd A = three_center_coulomb(basis, aux);
      const json& reference = record.at("three_center");
      REQUIRE(reference.size() == static_cast<std::size_t>(A.rows()));
      double three = 0.0;
      for (Eigen::Index r = 0; r < A.rows(); ++r)
        for (Eigen::Index c = 0; c < A.cols(); ++c)
          three = std::max(three, std::abs(A(r, c) - reference[static_cast<std::size_t>(r)]
                                                              [static_cast<std::size_t>(c)].get<double>()));
      INFO("two-centre " << two << ", three-centre " << three);
      CHECK(two <= 1e-10);
      CHECK(three <= 1e-10);
      CHECK((V - V.transpose()).cwiseAbs().maxCoeff() < 1e-13);
      SUCCEED(name << ": two-centre " << two << ", three-centre " << three);
    }
  }
}

namespace {

constexpr double acceptance_threshold = 1e-7;  // Hartree

struct Summary {
  int cases = 0;
  double energy = 0.0;
};

void compare_energies(bool slow_cases, Summary& summary) {
  const json reference = testing::load_reference("m3/df.json");
  const json& settings = reference.at("meta").at("settings");
  const BasisLibrary aux_library = BasisLibrary::load(settings.at("auxiliary_basis").get<std::string>());
  Method method;
  method.grid = {settings.at("grid").at("n_radial").get<int>(), settings.at("grid").at("n_angular").get<int>()};
  ScfOptions options;
  options.energy_tolerance = 1e-10;
  options.error_tolerance = 1e-7;

  for (const json& record : reference.at("results")) {
    const std::string label = record.at("molecule").get<std::string>() + " / " +
                              record.at("functional").get<std::string>();
    const Molecule molecule = testing::molecule_from_reference(record, label);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(record.at("basis").get<std::string>()),
                                           record.at("spherical").get<bool>());
    const bool quick = molecule.n_atoms() <= 2 && basis.n_primitives() <= 12;
    if (quick == slow_cases) continue;
    INFO(label);
    const BasisSet aux = BasisSet::build(molecule, aux_library, true);
    REQUIRE(aux.n_functions() == record.at("n_auxiliary").get<int>());

    DensityFittedJK jk(basis, aux);
    method.functional = Functional::by_name(record.at("functional").get<std::string>());
    const ScfResult scf = run_scf(molecule, basis, method, jk, options);
    CHECK(scf.converged);
    const double deviation = std::abs(scf.energy_total - record.at("energy_total").get<double>());
    CHECK(deviation <= acceptance_threshold);
    ++summary.cases;
    summary.energy = std::max(summary.energy, deviation);
    std::printf("DF %-8s %-4s %2d iterations  E = %.10f  |dE| = %.2e  DF error (reference) = %+.2e\n",
                record.at("molecule").get<std::string>().c_str(), record.at("functional").get<std::string>().c_str(),
                scf.iterations, scf.energy_total, deviation, record.at("fitting_error").get<double>());
    std::fflush(stdout);
  }
}

}  // namespace

TEST_CASE("M3: energies with density fitting equal the reference (atoms and H2)", "[m3][regression][df]") {
  Summary summary;
  compare_energies(false, summary);
  CHECK(summary.cases == 6);
  SUCCEED("largest deviation over " << summary.cases << " cases: " << summary.energy << " Hartree");
}

TEST_CASE("M3: energies with density fitting equal the reference (all other cases)",
          "[m3][regression][df][slow]") {
  Summary summary;
  compare_energies(true, summary);
  CHECK(summary.cases == 18);
  SUCCEED("largest deviation over " << summary.cases << " cases: " << summary.energy << " Hartree");
}
