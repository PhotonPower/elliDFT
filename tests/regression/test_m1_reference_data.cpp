// Milestone M0, third acceptance criterion: the reference data for M1 exist,
// are complete and internally consistent. These checks concern the stored
// reference files only (tests/reference/data/m1, produced by
// tests/reference/regenerate.sh m1); elliDFT computes none of the numbers yet.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <filesystem>
#include <set>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "reference_data.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using nlohmann::json;

namespace {

std::vector<std::string> reference_files(const std::string& prefix) {
  std::vector<std::string> names;
  const auto dir = testing::source_dir() / "tests" / "reference" / "data" / "m1";
  for (const auto& entry : std::filesystem::directory_iterator(dir)) {
    const std::string name = entry.path().filename().string();
    if (name.rfind(prefix, 0) == 0 && entry.path().extension() == ".json") names.push_back(name);
  }
  return names;
}

std::size_t n_pairs(std::size_t n) { return n * (n + 1) / 2; }

void check_meta(const json& reference) {
  const json& meta = reference.at("meta");
  CHECK(meta.at("program").get<std::string>() == "PySCF");
  CHECK_FALSE(meta.at("pyscf_version").get<std::string>().empty());
  CHECK(meta.at("bohr_in_angstrom").get<double>() == 0.529177210903);
  CHECK(meta.contains("settings"));
}

}  // namespace

TEST_CASE("M1 reference: Boys function table", "[m1][reference]") {
  const json reference = testing::load_reference("m1/boys.json");
  const auto x = reference.at("x").get<std::vector<double>>();
  const auto values = reference.at("values").get<std::vector<std::vector<double>>>();
  const std::size_t n_max = reference.at("meta").at("n_max").get<std::size_t>();

  REQUIRE(n_max >= 24);  // acceptance criterion M1: n <= 24
  REQUIRE(x.size() == values.size());
  REQUIRE(x.size() > 500);
  CHECK(x.front() == 0.0);
  CHECK(x.back() >= 500.0);  // acceptance criterion M1: 0 <= x <= 500

  for (std::size_t i = 0; i < x.size(); ++i) {
    INFO("x = " << x[i]);
    if (i > 0) CHECK(x[i] > x[i - 1]);
    REQUIRE(values[i].size() == n_max + 1);
    for (std::size_t n = 0; n <= n_max; ++n) {
      // 0 < F_n(x) <= 1/(2n+1), strictly decreasing in n
      CHECK(values[i][n] > 0.0);
      CHECK(values[i][n] <= 1.0 / static_cast<double>(2 * n + 1));
      if (n > 0) CHECK(values[i][n] < values[i][n - 1]);
    }
  }
  // F_n(0) = 1/(2n+1) exactly (correctly rounded)
  for (std::size_t n = 0; n <= n_max; ++n)
    CHECK(values[0][n] == 1.0 / static_cast<double>(2 * n + 1));
}

TEST_CASE("M1 reference: one-electron matrix files", "[m1][reference]") {
  const auto files = reference_files("ints1e_");
  REQUIRE(files.size() == 57);
  for (const std::string& name : files) {
    DYNAMIC_SECTION(name) {
      const json reference = testing::load_reference("m1/" + name);
      check_meta(reference);
      const std::size_t n = reference.at("n_functions").get<std::size_t>();
      CHECK(reference.at("function_labels").size() == n);
      CHECK(reference.at("molecule").at("symbols").size() ==
            reference.at("molecule").at("geometry_bohr").size());
      for (const char* key : {"overlap", "kinetic", "nuclear", "dipole_x", "dipole_y", "dipole_z"}) {
        INFO(key);
        CHECK(reference.at("matrices").at(key).size() == n_pairs(n));
      }
    }
  }
}

TEST_CASE("M1 reference: two-electron integral files", "[m1][reference]") {
  const auto files = reference_files("eri_");
  REQUIRE(files.size() == 11);
  for (const std::string& name : files) {
    DYNAMIC_SECTION(name) {
      const json reference = testing::load_reference("m1/" + name);
      check_meta(reference);
      const std::size_t n = reference.at("n_functions").get<std::size_t>();
      CHECK(reference.at("eri_packed").size() == n_pairs(n_pairs(n)));
    }
  }
}

TEST_CASE("M1 reference: high angular momentum cases", "[m1][reference]") {
  for (const std::string name : {"highl_sph.json", "highl_cart.json"}) {
    DYNAMIC_SECTION(name) {
      const json reference = testing::load_reference("m1/" + name);
      check_meta(reference);
      const std::size_t n = reference.at("n_functions").get<std::size_t>();
      CHECK(reference.at("matrices").at("overlap").size() == n_pairs(n));

      // Angular momenta 0..6 are all present.
      std::set<int> ls;
      std::size_t n_shells = 0;
      for (const json& atom : reference.at("basis").at("shells_per_atom"))
        for (const json& shell : atom) {
          ls.insert(shell.at("l").get<int>());
          ++n_shells;
        }
      CHECK(ls == std::set<int>{0, 1, 2, 3, 4, 5, 6});
      CHECK(n_shells == reference.at("n_shells").get<std::size_t>());

      REQUIRE_FALSE(reference.at("eri_shell_quartets").empty());
      for (const json& quartet : reference.at("eri_shell_quartets")) {
        const auto shape = quartet.at("shape").get<std::vector<std::size_t>>();
        REQUIRE(shape.size() == 4);
        CHECK(quartet.at("values").size() == shape[0] * shape[1] * shape[2] * shape[3]);
        for (const json& shell : quartet.at("shells")) CHECK(shell.get<std::size_t>() < n_shells);
      }
    }
  }
}

TEST_CASE("M1 reference: Hartree-Fock results", "[m1][reference]") {
  const json reference = testing::load_reference("m1/hf.json");
  check_meta(reference);
  // Reference SCF convergence as demanded in the specification, section 7.
  CHECK(reference.at("meta").at("settings").at("conv_tol").get<double>() <= 1e-12);

  const json& results = reference.at("results");
  // 12 molecules (H2 ... benzene) in STO-3G, 6-31G*, def2-SVP, cc-pVDZ, aug-cc-pVDZ
  REQUIRE(results.size() == 60);

  std::set<std::string> molecules, basis_sets;
  for (const json& r : results) {
    INFO(r.at("molecule").get<std::string>() << " / " << r.at("basis").get<std::string>());
    molecules.insert(r.at("molecule").get<std::string>());
    basis_sets.insert(r.at("basis").get<std::string>());

    const double total = r.at("energy_total").get<double>();
    const double parts = r.at("energy_nuclear_repulsion").get<double>() +
                         r.at("energy_one_electron").get<double>() +
                         r.at("energy_two_electron").get<double>();
    CHECK_THAT(total - parts, WithinAbs(0.0, 1e-10));
    CHECK(total < 0.0);

    // Closed shell, no functions removed, Tr(DS) = number of electrons.
    const int n_electrons = r.at("n_electrons").get<int>();
    CHECK(n_electrons % 2 == 0);
    CHECK(r.at("n_orbitals").get<int>() == r.at("n_functions").get<int>());
    CHECK(r.at("orbital_energies").size() == r.at("n_functions").get<std::size_t>());
    CHECK_THAT(r.at("trace_DS").get<double>(), WithinAbs(n_electrons, 1e-9));
    CHECK(r.at("overlap_min_eigenvalue").get<double>() > 0.0);
    CHECK(r.at("dipole_au").size() == 3);
    // 6-31G* is cartesian, the other basis sets are spherical.
    CHECK(r.at("spherical").get<bool>() == (r.at("basis").get<std::string>() != "6-31G*"));
  }
  CHECK(molecules.size() == 12);
  CHECK(basis_sets ==
        std::set<std::string>{"STO-3G", "6-31G*", "def2-SVP", "cc-pVDZ", "aug-cc-pVDZ"});
}

// External anchor named in the specification (M1, attributed there to
// Szabo/Ostlund): H2 in STO-3G at 1.4 Bohr gives -1.1167 Hartree. Here the
// anchor checks the reference data; in M1 it checks elliDFT itself.
TEST_CASE("M1 reference: H2 anchor of the specification", "[m1][reference]") {
  const json reference = testing::load_reference("m1/hf.json");
  bool found = false;
  for (const json& r : reference.at("results")) {
    if (r.at("molecule") != "h2" || r.at("basis") != "STO-3G") continue;
    found = true;
    const auto geometry = r.at("geometry_bohr").get<std::vector<std::vector<double>>>();
    REQUIRE(geometry.size() == 2);
    const double distance = std::hypot(geometry[1][0] - geometry[0][0], geometry[1][1] - geometry[0][1],
                                       geometry[1][2] - geometry[0][2]);
    CHECK_THAT(distance, WithinAbs(1.4, 1e-14));
    // -1.1167 is given to four decimals.
    CHECK_THAT(r.at("energy_total").get<double>(), WithinAbs(-1.1167, 5e-5));
  }
  CHECK(found);
}
