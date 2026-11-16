#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <stdexcept>

#include "core/elements.hpp"
#include "core/molecule.hpp"
#include "core/units.hpp"

using namespace ellidft;
using Catch::Matchers::WithinRel;

TEST_CASE("unit conversions are mutually inverse", "[core][units]") {
  // Values fixed by the specification (CODATA 2018).
  CHECK(units::bohr_in_angstrom == 0.529177210903);
  CHECK(units::hartree_in_ev == 27.211386245988);
  CHECK_THAT(units::angstrom_to_bohr(units::bohr_in_angstrom), WithinRel(1.0, 1e-15));
  CHECK_THAT(units::bohr_to_angstrom(units::angstrom_to_bohr(1.2345)), WithinRel(1.2345, 1e-15));
}

TEST_CASE("element symbols and atomic numbers", "[core][elements]") {
  CHECK(atomic_number("H") == 1);
  CHECK(atomic_number("He") == 2);
  CHECK(atomic_number("C") == 6);
  CHECK(atomic_number("Cl") == 17);
  CHECK(atomic_number("Ar") == 18);

  SECTION("matching ignores case") {
    CHECK(atomic_number("cl") == 17);
    CHECK(atomic_number("CL") == 17);
    CHECK(atomic_number("h") == 1);
  }
  SECTION("unknown and unsupported symbols give 0") {
    CHECK(atomic_number("") == 0);
    CHECK(atomic_number("X") == 0);
    CHECK(atomic_number("K") == 0);   // beyond Ar
    CHECK(atomic_number("Br") == 0);
    CHECK(atomic_number("H1") == 0);
  }
  SECTION("round trip for all supported elements") {
    for (int z = 1; z <= max_atomic_number; ++z) CHECK(atomic_number(element_symbol(z)) == z);
  }
  SECTION("out of range atomic numbers throw") {
    CHECK_THROWS_AS(element_symbol(0), std::out_of_range);
    CHECK_THROWS_AS(element_symbol(max_atomic_number + 1), std::out_of_range);
  }
}

TEST_CASE("molecule electron count", "[core][molecule]") {
  Molecule water;
  water.atoms = {{8, {0.0, 0.0, 0.0}}, {1, {0.0, 1.0, 1.0}}, {1, {0.0, -1.0, 1.0}}};
  CHECK(water.n_atoms() == 3);
  CHECK(water.nuclear_charge() == 10);
  CHECK(water.n_electrons() == 10);

  water.charge = 1;
  CHECK(water.n_electrons() == 9);
  water.charge = -1;
  CHECK(water.n_electrons() == 11);
}
