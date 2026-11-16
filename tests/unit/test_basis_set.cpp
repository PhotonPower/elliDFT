#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <sstream>
#include <string>

#include "basis/basis_set.hpp"
#include "core/molecule.hpp"

using namespace ellidft;
using Catch::Matchers::ContainsSubstring;

namespace {

// Small artificial basis: H has s, p; O has s, s, p, d, f.
const char* const test_basis = R"(
BASIS "ao basis" SPHERICAL
H    S
      3.0   0.2
      0.5   0.8
H    P
      1.0   1.0
O    SP
      8.0   0.1   0.3
      1.5   0.9   0.7
O    S
      0.3   1.0
O    D
      1.2   1.0
O    F
      1.4   1.0
END
)";

BasisLibrary test_library() {
  std::istringstream in(test_basis);
  return BasisLibrary::parse_nwchem(in, "test");
}

Molecule water() {
  Molecule mol;
  mol.name = "water";
  mol.atoms = {{8, {0.0, 0.0, 0.2}}, {1, {0.0, 1.4, -0.9}}, {1, {0.0, -1.4, -0.9}}};
  return mol;
}

}  // namespace

TEST_CASE("number of functions per shell", "[basis]") {
  const int cartesian[] = {1, 3, 6, 10, 15, 21, 28};
  const int spherical[] = {1, 3, 5, 7, 9, 11, 13};
  for (int l = 0; l <= max_angular_momentum; ++l) {
    CHECK(n_cartesian(l) == cartesian[l]);
    CHECK(n_spherical(l) == spherical[l]);
    // Direct count of the monomials x^i y^j z^k with i + j + k = l.
    int count = 0;
    for (int i = 0; i <= l; ++i)
      for (int j = 0; j <= l - i; ++j) ++count;
    CHECK(n_cartesian(l) == count);
  }
}

TEST_CASE("basis of a molecule: spherical", "[basis]") {
  const Molecule mol = water();
  const BasisSet basis = BasisSet::build(mol, test_library());

  CHECK(basis.name() == "test");
  CHECK(basis.spherical());
  // O: s s p d f (5 shells), H: s p (2 shells each)
  REQUIRE(basis.n_shells() == 9);
  CHECK(basis.max_l() == 3);
  // O: 1 + 1 + 3 + 5 + 7 = 17, H: 1 + 3 = 4
  CHECK(basis.n_functions() == 17 + 4 + 4);
  // O: 2 + 1 + 2 + 1 + 1 = 7 primitives, H: 2 + 1 = 3
  CHECK(basis.n_primitives() == 7 + 3 + 3);

  const int expected_atom[] = {0, 0, 0, 0, 0, 1, 1, 2, 2};
  const int expected_l[] = {0, 0, 1, 2, 3, 0, 1, 0, 1};
  const int expected_first[] = {0, 1, 2, 5, 10, 17, 18, 21, 22};
  for (std::size_t i = 0; i < basis.n_shells(); ++i) {
    INFO("shell " << i);
    const Shell& shell = basis.shells()[i];
    CHECK(shell.atom == expected_atom[i]);
    CHECK(shell.l == expected_l[i]);
    CHECK(shell.spherical);
    CHECK(shell.center == mol.atoms[static_cast<std::size_t>(shell.atom)].r);
    CHECK(basis.first_function(i) == expected_first[i]);
  }
  CHECK_THROWS(basis.first_function(basis.n_shells()));

  // Data of the SP block arrives in both shells.
  CHECK(basis.shells()[0].exponents == std::vector<double>{8.0, 1.5});
  CHECK(basis.shells()[0].coefficients == std::vector<double>{0.1, 0.9});
  CHECK(basis.shells()[2].exponents == std::vector<double>{8.0, 1.5});
  CHECK(basis.shells()[2].coefficients == std::vector<double>{0.3, 0.7});
}

TEST_CASE("basis of a molecule: cartesian override", "[basis]") {
  const BasisSet basis = BasisSet::build(water(), test_library(), /*spherical=*/false);
  CHECK_FALSE(basis.spherical());
  REQUIRE(basis.n_shells() == 9);
  // O: 1 + 1 + 3 + 6 + 10 = 21, H: 4
  CHECK(basis.n_functions() == 21 + 4 + 4);
  for (const Shell& shell : basis.shells()) CHECK_FALSE(shell.spherical);
  CHECK(basis.first_function(4) == 11);
  CHECK(basis.first_function(5) == 21);

  // Explicit choice equal to the file default changes nothing.
  CHECK(BasisSet::build(water(), test_library(), true).n_functions() == 25);
}

TEST_CASE("basis of a molecule: missing element", "[basis]") {
  Molecule mol = water();
  mol.atoms.push_back({6, {5.0, 0.0, 0.0}});
  CHECK_THROWS_WITH(BasisSet::build(mol, test_library()),
                    ContainsSubstring("no functions for C (atom 4)"));
}

TEST_CASE("basis of an empty molecule is empty", "[basis]") {
  const BasisSet basis = BasisSet::build(Molecule{}, test_library());
  CHECK(basis.n_shells() == 0);
  CHECK(basis.n_functions() == 0);
  CHECK(basis.n_primitives() == 0);
}
