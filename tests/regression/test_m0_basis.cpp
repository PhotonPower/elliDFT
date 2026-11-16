// Milestone M0: basis set structure against the reference code.
//
// Reference: tests/reference/data/m0/*.json (PySCF, tests/reference/gen_m0.py).

#include <catch2/catch_test_macros.hpp>

#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/elements.hpp"
#include "core/xyz.hpp"
#include "reference_data.hpp"

using namespace ellidft;
using nlohmann::json;

// Exponents and coefficients are read from the same files by both programs
// and must agree exactly: decimal-to-double conversion is correctly rounded on
// both sides, so no tolerance is involved.
TEST_CASE("M0: shells of all elements and basis sets equal the reference", "[m0][regression]") {
  const json reference = testing::load_reference("m0/basis_shells.json");
  const json& basis_sets = reference.at("basis_sets");
  REQUIRE(basis_sets.size() == 11);

  for (const auto& [basis_name, record] : basis_sets.items()) {
    DYNAMIC_SECTION(basis_name) {
      const BasisLibrary library = BasisLibrary::load(basis_name);
      CHECK(library.spherical() == record.at("spherical").get<bool>());
      REQUIRE(record.at("elements").size() == static_cast<std::size_t>(max_atomic_number));

      for (const auto& [symbol, ref_shells] : record.at("elements").items()) {
        INFO("element " << symbol);
        const int Z = atomic_number(symbol);
        REQUIRE(Z > 0);
        REQUIRE(library.has_element(Z));
        const auto& shells = library.shells(Z);
        REQUIRE(shells.size() == ref_shells.size());
        for (std::size_t i = 0; i < shells.size(); ++i) {
          INFO("shell " << i);
          CHECK(shells[i].l == ref_shells[i].at("l").get<int>());
          CHECK(shells[i].exponents == ref_shells[i].at("exponents").get<std::vector<double>>());
          CHECK(shells[i].coefficients ==
                ref_shells[i].at("coefficients").get<std::vector<double>>());
        }
      }
    }
  }
}

// Acceptance criterion M0: number of shells and basis functions for all test
// molecules, all basis sets, spherical and cartesian.
TEST_CASE("M0: basis sizes of all test molecules equal the reference", "[m0][regression]") {
  const json reference = testing::load_reference("m0/molecules.json");
  const json& molecules = reference.at("molecules");
  // 12 small molecules, 2 repeat units, 7 styrene and 3 MMA oligomers, and the
  // isotactic control series (2 x 3 molecules)
  REQUIRE(molecules.size() == 30);

  std::size_t n_compared = 0;
  for (const auto& [name, record] : molecules.items()) {
    DYNAMIC_SECTION(name) {
      const Molecule mol = testing::molecule_from_reference(record, name);
      CHECK(mol.n_electrons() == record.at("n_electrons").get<int>());
      CHECK(mol.n_electrons() % 2 == 0);  // closed shell

      REQUIRE(record.at("basis").size() == 11);
      for (const auto& [basis_name, sizes] : record.at("basis").items()) {
        const BasisLibrary library = BasisLibrary::load(basis_name);
        for (const bool spherical : {true, false}) {
          INFO(basis_name << (spherical ? " spherical" : " cartesian"));
          const json& ref = sizes.at(spherical ? "spherical" : "cartesian");
          const BasisSet basis = BasisSet::build(mol, library, spherical);

          CHECK(basis.n_functions() == ref.at("n_functions").get<int>());
          CHECK(basis.n_shells() == ref.at("n_shells").get<std::size_t>());
          CHECK(basis.n_primitives() == ref.at("n_primitives").get<std::size_t>());

          std::vector<int> per_l(static_cast<std::size_t>(basis.max_l() + 1), 0);
          for (const Shell& shell : basis.shells()) ++per_l[static_cast<std::size_t>(shell.l)];
          CHECK(per_l == ref.at("n_shells_per_l").get<std::vector<int>>());
          ++n_compared;
        }
      }
    }
  }
  // Catch2 re-runs the test body once per dynamic section, so only the
  // comparisons of one molecule are counted in each pass.
  CHECK(n_compared == 22);
}

// The XYZ reader must reproduce the Bohr coordinates that the reference
// scripts pass to PySCF (same constant, same operation: bit-identical).
TEST_CASE("M0: XYZ files give the reference coordinates in Bohr", "[m0][regression]") {
  const json reference = testing::load_reference("m0/molecules.json");
  for (const auto& [name, record] : reference.at("molecules").items()) {
    DYNAMIC_SECTION(name) {
      const Molecule from_reference = testing::molecule_from_reference(record, name);
      const Molecule from_xyz =
          read_xyz_file(testing::source_dir() / record.at("xyz_file").get<std::string>());
      REQUIRE(from_xyz.n_atoms() == from_reference.n_atoms());
      CHECK(from_xyz.name == record.at("comment").get<std::string>());
      for (std::size_t i = 0; i < from_xyz.n_atoms(); ++i) {
        INFO("atom " << i);
        CHECK(from_xyz.atoms[i].Z == from_reference.atoms[i].Z);
        CHECK(from_xyz.atoms[i].r == from_reference.atoms[i].r);
      }
    }
  }
}
