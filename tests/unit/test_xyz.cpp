#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <sstream>
#include <stdexcept>
#include <string>

#include "core/units.hpp"
#include "core/xyz.hpp"

using namespace ellidft;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Molecule parse(const std::string& text, LengthUnit unit = LengthUnit::Angstrom) {
  std::istringstream in(text);
  return read_xyz(in, unit);
}

}  // namespace

TEST_CASE("XYZ reader converts Angstrom to Bohr", "[core][xyz]") {
  const Molecule mol = parse("3\n  water, test  \n"
                             "O   0.000000   0.000000   0.117790\n"
                             "H   0.000000   0.755453  -0.471161\n"
                             "h   0.000000  -0.755453  -0.471161\n");
  REQUIRE(mol.n_atoms() == 3);
  CHECK(mol.name == "water, test");
  CHECK(mol.charge == 0);
  CHECK(mol.atoms[0].Z == 8);
  CHECK(mol.atoms[1].Z == 1);
  CHECK(mol.atoms[2].Z == 1);
  CHECK(mol.atoms[0].r[0] == 0.0);
  CHECK_THAT(mol.atoms[0].r[2], WithinRel(0.117790 / units::bohr_in_angstrom, 1e-15));
  CHECK_THAT(mol.atoms[1].r[1], WithinRel(0.755453 / units::bohr_in_angstrom, 1e-15));
  CHECK_THAT(mol.atoms[2].r[1], WithinRel(-0.755453 / units::bohr_in_angstrom, 1e-15));
  CHECK_THAT(mol.atoms[2].r[2], WithinRel(-0.471161 / units::bohr_in_angstrom, 1e-15));
}

TEST_CASE("XYZ reader takes Bohr unchanged", "[core][xyz]") {
  const Molecule mol = parse("2\nH2\nH 0 0 0\nH 0 0 1.4\n", LengthUnit::Bohr);
  REQUIRE(mol.n_atoms() == 2);
  CHECK(mol.atoms[1].r[2] == 1.4);
  CHECK(mol.n_electrons() == 2);
}

TEST_CASE("XYZ reader tolerates harmless variations", "[core][xyz]") {
  SECTION("empty comment line, exponent notation, extra columns, trailing blank lines") {
    const Molecule mol = parse("  2  \n\nNe 1.0e0 -2.5E-1 +3 extra columns\nAr\t0\t0\t0\n\n  \n");
    REQUIRE(mol.n_atoms() == 2);
    CHECK(mol.name.empty());
    CHECK(mol.atoms[0].Z == 10);
    CHECK(mol.atoms[1].Z == 18);
    CHECK_THAT(mol.atoms[0].r[0], WithinRel(1.0 / units::bohr_in_angstrom, 1e-15));
    CHECK_THAT(mol.atoms[0].r[1], WithinRel(-0.25 / units::bohr_in_angstrom, 1e-15));
    CHECK_THAT(mol.atoms[0].r[2], WithinRel(3.0 / units::bohr_in_angstrom, 1e-15));
    CHECK_THAT(mol.atoms[1].r[0], WithinAbs(0.0, 0.0));
  }
  SECTION("Windows line endings") {
    const Molecule mol = parse("1\r\nhelium\r\nHe 0.0 0.0 0.0\r\n");
    REQUIRE(mol.n_atoms() == 1);
    CHECK(mol.name == "helium");
    CHECK(mol.atoms[0].Z == 2);
  }
}

TEST_CASE("XYZ reader rejects malformed input", "[core][xyz]") {
  CHECK_THROWS_WITH(parse(""), ContainsSubstring("atom count"));
  CHECK_THROWS_WITH(parse("abc\ncomment\nH 0 0 0\n"), ContainsSubstring("atom count"));
  CHECK_THROWS_WITH(parse("0\ncomment\n"), ContainsSubstring("atom count"));
  CHECK_THROWS_WITH(parse("-1\ncomment\n"), ContainsSubstring("atom count"));
  CHECK_THROWS_WITH(parse("1.5\ncomment\nH 0 0 0\n"), ContainsSubstring("atom count"));
  CHECK_THROWS_WITH(parse("2\n"), ContainsSubstring("comment"));

  SECTION("too few atoms") {
    CHECK_THROWS_WITH(parse("2\ncomment\nH 0 0 0\n"), ContainsSubstring("file ends after 1 of 2"));
  }
  SECTION("too many atoms") {
    CHECK_THROWS_WITH(parse("1\ncomment\nH 0 0 0\nH 0 0 1\n"),
                      ContainsSubstring("unexpected content"));
  }
  SECTION("unknown and unsupported elements") {
    CHECK_THROWS_WITH(parse("1\ncomment\nXx 0 0 0\n"), ContainsSubstring("element 'Xx'"));
    CHECK_THROWS_WITH(parse("1\ncomment\nBr 0 0 0\n"), ContainsSubstring("element 'Br'"));
    CHECK_THROWS_WITH(parse("1\ncomment\n6 0 0 0\n"), ContainsSubstring("element '6'"));
  }
  SECTION("bad coordinates") {
    CHECK_THROWS_WITH(parse("1\ncomment\nH 0 0\n"), ContainsSubstring("<symbol> <x> <y> <z>"));
    CHECK_THROWS_WITH(parse("1\ncomment\nH 0 zero 0\n"), ContainsSubstring("coordinate 'zero'"));
    CHECK_THROWS_WITH(parse("1\ncomment\nH 0 1.0x 0\n"), ContainsSubstring("coordinate '1.0x'"));
    CHECK_THROWS_WITH(parse("1\ncomment\nH 0 nan 0\n"), ContainsSubstring("coordinate 'nan'"));
    CHECK_THROWS_WITH(parse("1\ncomment\nH 0 inf 0\n"), ContainsSubstring("coordinate 'inf'"));
  }
  SECTION("error message carries the line number") {
    CHECK_THROWS_WITH(parse("2\ncomment\nH 0 0 0\nH 0 0 x\n"), ContainsSubstring("line 4"));
  }
}

TEST_CASE("XYZ file reader reports missing files", "[core][xyz]") {
  CHECK_THROWS_WITH(read_xyz_file("/nonexistent/molecule.xyz"), ContainsSubstring("cannot open"));
}
