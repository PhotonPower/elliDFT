#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "basis/basis_library.hpp"
#include "core/elements.hpp"

using namespace ellidft;
using Catch::Matchers::ContainsSubstring;

namespace {

BasisLibrary parse(const std::string& text) {
  std::istringstream in(text);
  return BasisLibrary::parse_nwchem(in, "test");
}

// All basis sets fixed in the specification, section 2.
const std::vector<std::string> all_basis_sets = {
    "STO-3G",     "6-31G*",  "def2-SVP",    "def2-TZVP",   "def2-SVPD",            "def2-TZVPD",
    "cc-pVDZ",    "aug-cc-pVDZ", "aug-cc-pVTZ", "def2-universal-JKFIT", "def2-universal-JFIT"};

using ShellCounts = std::array<int, max_angular_momentum + 1>;

// Parses the contracted pattern of a Basis Set Exchange comment line
//   "#BASIS SET: (9s,4p,1d) -> [3s,2p,1d]"
// into the number of contracted shells per angular momentum.
ShellCounts parse_bse_pattern(const std::string& line) {
  ShellCounts counts{};
  const auto open = line.find('[');
  const auto close = line.find(']');
  REQUIRE(open != std::string::npos);
  REQUIRE(close != std::string::npos);
  std::string body = line.substr(open + 1, close - open - 1);
  std::replace(body.begin(), body.end(), ',', ' ');
  std::istringstream in(body);
  std::string item;
  while (in >> item) {
    const int l = angular_momentum_from_letter(item.back());
    REQUIRE(l >= 0);
    counts[static_cast<std::size_t>(l)] += std::stoi(item.substr(0, item.size() - 1));
  }
  return counts;
}

// Contracted shell pattern per element, read from the comments that the Basis
// Set Exchange writes in front of each element. Independent of the parser
// under test.
std::map<int, ShellCounts> bse_patterns(const std::string& basis_name) {
  std::ifstream in(basis_directory() / basis_file_name(basis_name));
  REQUIRE(in.good());
  std::map<int, ShellCounts> out;
  std::string line;
  bool have_pattern = false;
  ShellCounts pattern{};
  while (std::getline(in, line)) {
    if (line.rfind("#BASIS SET:", 0) == 0) {
      pattern = parse_bse_pattern(line);
      have_pattern = true;
    } else if (have_pattern && !line.empty() &&
               std::isalpha(static_cast<unsigned char>(line[0]))) {
      std::istringstream header(line);
      std::string symbol;
      header >> symbol;
      const int Z = atomic_number(symbol);
      if (Z > 0) out[Z] = pattern;
      have_pattern = false;
    }
  }
  return out;
}

}  // namespace

TEST_CASE("angular momentum letters", "[basis]") {
  CHECK(angular_momentum_from_letter('s') == 0);
  CHECK(angular_momentum_from_letter('P') == 1);
  CHECK(angular_momentum_from_letter('d') == 2);
  CHECK(angular_momentum_from_letter('F') == 3);
  CHECK(angular_momentum_from_letter('g') == 4);
  CHECK(angular_momentum_from_letter('h') == 5);
  CHECK(angular_momentum_from_letter('I') == 6);
  CHECK(angular_momentum_from_letter('x') == -1);
  CHECK(angular_momentum_from_letter('l') == -1);
  for (int l = 0; l <= max_angular_momentum; ++l)
    CHECK(angular_momentum_from_letter(angular_momentum_letter(l)) == l);
  CHECK_THROWS_AS(angular_momentum_letter(-1), std::out_of_range);
  CHECK_THROWS_AS(angular_momentum_letter(max_angular_momentum + 1), std::out_of_range);
}

TEST_CASE("basis file names", "[basis]") {
  CHECK(basis_file_name("STO-3G") == "sto-3g.nw");
  CHECK(basis_file_name("6-31G*") == "6-31g_st.nw");
  CHECK(basis_file_name("def2-universal-JKFIT") == "def2-universal-jkfit.nw");
}

TEST_CASE("NWChem parser: segmented shells", "[basis][parser]") {
  const BasisLibrary lib = parse(R"(
# a comment
BASIS "ao basis" SPHERICAL PRINT
#BASIS SET: (3s) -> [1s]
H    S
      0.3425250914E+01       0.1543289673E+00
      0.6239137298E+00       0.5353281423E+00
      0.1688554040E+00       0.4446345422E+00
H    P
      1.1D+00   1.0
END
)");
  CHECK(lib.name() == "test");
  CHECK(lib.spherical());
  REQUIRE(lib.has_element(1));
  CHECK_FALSE(lib.has_element(2));
  CHECK_THROWS_AS(lib.shells(2), std::out_of_range);

  const auto& shells = lib.shells(1);
  REQUIRE(shells.size() == 2);
  CHECK(shells[0].l == 0);
  CHECK(shells[0].exponents == std::vector<double>{3.425250914, 0.6239137298, 0.1688554040});
  CHECK(shells[0].coefficients == std::vector<double>{0.1543289673, 0.5353281423, 0.4446345422});
  CHECK(shells[1].l == 1);
  // Fortran exponent marker
  CHECK(shells[1].exponents == std::vector<double>{1.1});
  CHECK(shells[1].coefficients == std::vector<double>{1.0});
}

TEST_CASE("NWChem parser: combined SP shells are split", "[basis][parser]") {
  const BasisLibrary lib = parse(R"(
BASIS "ao basis" CARTESIAN PRINT
C    S
      10.0   1.0
C    SP
      2.0   -0.1   0.2
      0.5    0.6   0.7
C    D
      0.8    1.0
END
)");
  CHECK_FALSE(lib.spherical());
  const auto& shells = lib.shells(6);
  REQUIRE(shells.size() == 4);
  // ordered by angular momentum: s, s (from SP), p (from SP), d
  CHECK(shells[0].l == 0);
  CHECK(shells[0].exponents == std::vector<double>{10.0});
  CHECK(shells[1].l == 0);
  CHECK(shells[1].exponents == std::vector<double>{2.0, 0.5});
  CHECK(shells[1].coefficients == std::vector<double>{-0.1, 0.6});
  CHECK(shells[2].l == 1);
  CHECK(shells[2].exponents == std::vector<double>{2.0, 0.5});
  CHECK(shells[2].coefficients == std::vector<double>{0.2, 0.7});
  CHECK(shells[3].l == 2);
}

TEST_CASE("NWChem parser: general contractions become one shell per column", "[basis][parser]") {
  const BasisLibrary lib = parse(R"(
BASIS "ao basis" SPHERICAL
O    P
      5.0    0.1    0.0
      1.0    0.5    0.0
      0.2    0.6    1.0
O    S
      9.0    0.3   -0.1    0.0
      2.0    0.7   -0.2    0.0
      0.4    0.0    1.1    1.0
END
)");
  const auto& shells = lib.shells(8);
  REQUIRE(shells.size() == 5);
  // s shells first although the file lists p first; column order kept;
  // primitives with zero coefficient dropped.
  CHECK(shells[0].l == 0);
  CHECK(shells[0].exponents == std::vector<double>{9.0, 2.0});
  CHECK(shells[0].coefficients == std::vector<double>{0.3, 0.7});
  CHECK(shells[1].l == 0);
  CHECK(shells[1].exponents == std::vector<double>{9.0, 2.0, 0.4});
  CHECK(shells[1].coefficients == std::vector<double>{-0.1, -0.2, 1.1});
  CHECK(shells[2].l == 0);
  CHECK(shells[2].exponents == std::vector<double>{0.4});
  CHECK(shells[2].coefficients == std::vector<double>{1.0});
  CHECK(shells[3].l == 1);
  CHECK(shells[3].exponents == std::vector<double>{5.0, 1.0, 0.2});
  CHECK(shells[4].l == 1);
  CHECK(shells[4].exponents == std::vector<double>{0.2});
}

TEST_CASE("NWChem parser: header keywords and unsupported elements", "[basis][parser]") {
  SECTION("no keyword means cartesian (NWChem default)") {
    CHECK_FALSE(parse("BASIS \"ao basis\"\nH S\n 1.0 1.0\nEND\n").spherical());
  }
  SECTION("keywords are case-insensitive") {
    CHECK(parse("basis \"ao basis\" spherical\nh s\n 1.0 1.0\nend\n").spherical());
  }
  SECTION("elements beyond Ar are skipped") {
    const BasisLibrary lib =
        parse("BASIS \"ao basis\" SPHERICAL\nH S\n 1.0 1.0\nK S\n 2.0 1.0\n 1.0 0.5\nEND\n");
    CHECK(lib.has_element(1));
    CHECK(lib.shells(1).size() == 1);
    for (int z = 2; z <= max_atomic_number; ++z) CHECK_FALSE(lib.has_element(z));
  }
  SECTION("high angular momentum up to i") {
    const BasisLibrary lib = parse("BASIS \"x\" SPHERICAL\nH I\n 1.0 1.0\nH H\n 1.0 1.0\nEND\n");
    REQUIRE(lib.shells(1).size() == 2);
    CHECK(lib.shells(1)[0].l == 5);
    CHECK(lib.shells(1)[1].l == 6);
  }
}

TEST_CASE("NWChem parser rejects malformed input", "[basis][parser]") {
  CHECK_THROWS_WITH(parse(""), ContainsSubstring("no BASIS block"));
  CHECK_THROWS_WITH(parse("H S\n 1.0 1.0\n"), ContainsSubstring("outside of a BASIS block"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 1.0 1.0\n"), ContainsSubstring("missing END"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 1.0 1.0\nEND\nBASIS \"b\"\nEND\n"),
                    ContainsSubstring("more than one BASIS block"));
  CHECK_THROWS_WITH(parse("ECP\n H nelec 2\nEND\n"), ContainsSubstring("effective core potentials"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\n 1.0 1.0\nEND\n"), ContainsSubstring("before any shell"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\nEND\n"), ContainsSubstring("without primitives"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\nH P\n 1.0 1.0\nEND\n"),
                    ContainsSubstring("without primitives"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH X\n 1.0 1.0\nEND\n"), ContainsSubstring("shell type 'X'"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH\n 1.0 1.0\nEND\n"), ContainsSubstring("<element> <shell type>"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 1.0\nEND\n"), ContainsSubstring("at least one coefficient"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 1.0 abc\nEND\n"), ContainsSubstring("invalid number"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n -1.0 1.0\nEND\n"), ContainsSubstring("positive"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 0.0 1.0\nEND\n"), ContainsSubstring("positive"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 2.0 1.0 0.5\n 1.0 1.0\nEND\n"),
                    ContainsSubstring("inconsistent number"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH SP\n 2.0 1.0\nEND\n"),
                    ContainsSubstring("one coefficient column per letter"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 2.0 1.0 0.0\n 1.0 0.5 0.0\nEND\n"),
                    ContainsSubstring("only zero coefficients"));
  CHECK_THROWS_WITH(parse("BASIS \"a\"\nH S\n 1.0 1.0\n\n 2.0 x\nEND\n"), ContainsSubstring("line 5"));
  CHECK_THROWS_WITH(BasisLibrary::load("no-such-basis"), ContainsSubstring("cannot open"));
}

// Acceptance criterion M0: all basis set files of section 2 load for H to Cl
// (Ar in addition, see core/elements.hpp).
TEST_CASE("all basis set files load for H to Ar", "[basis][files]") {
  for (const std::string& name : all_basis_sets) {
    DYNAMIC_SECTION(name) {
      const BasisLibrary lib = BasisLibrary::load(name);
      CHECK(lib.name() == name);
      // 6-31G* is defined with cartesian d functions, all others spherical.
      CHECK(lib.spherical() == (name != "6-31G*"));

      for (int Z = 1; Z <= max_atomic_number; ++Z) {
        INFO("element " << element_symbol(Z));
        REQUIRE(lib.has_element(Z));
        const auto& shells = lib.shells(Z);
        REQUIRE_FALSE(shells.empty());
        CHECK(shells.front().l == 0);
        int previous_l = 0;
        for (const ShellDefinition& shell : shells) {
          CHECK(shell.l >= previous_l);  // ordered by angular momentum
          CHECK(shell.l <= max_angular_momentum);
          previous_l = shell.l;
          REQUIRE(shell.n_primitives() >= 1);
          REQUIRE(shell.coefficients.size() == shell.exponents.size());
          for (std::size_t k = 0; k < shell.n_primitives(); ++k) {
            CHECK(shell.exponents[k] > 0.0);
            CHECK(shell.coefficients[k] != 0.0);
          }
        }
      }
    }
  }
}

// Check without the reference code: the number of contracted shells per
// angular momentum must equal the pattern "[3s,2p,1d]" that the Basis Set
// Exchange writes as a comment in front of each element.
TEST_CASE("shell counts agree with the Basis Set Exchange patterns", "[basis][files]") {
  for (const std::string& name : all_basis_sets) {
    DYNAMIC_SECTION(name) {
      const BasisLibrary lib = BasisLibrary::load(name);
      const auto patterns = bse_patterns(name);
      REQUIRE(patterns.size() == static_cast<std::size_t>(max_atomic_number));
      for (const auto& [Z, expected] : patterns) {
        INFO("element " << element_symbol(Z));
        ShellCounts counts{};
        for (const ShellDefinition& shell : lib.shells(Z)) ++counts[static_cast<std::size_t>(shell.l)];
        CHECK(counts == expected);
      }
    }
  }
}
