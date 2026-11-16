// Two-electron integrals: checks that need no reference code.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <sstream>
#include <vector>

#include "basis/basis_set.hpp"
#include "full_eri.hpp"
#include "integrals/shell_pair.hpp"
#include "integrals/two_electron.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Shell s_shell(int atom, const Vec3& center, double exponent) {
  Shell shell;
  shell.atom = atom;
  shell.center = center;
  shell.l = 0;
  shell.exponents = {exponent};
  shell.coefficients = {1.0};
  return shell;
}

const char* const test_basis = R"(
BASIS "ao basis" SPHERICAL
H    S
      2.9   0.3
      0.5   0.8
H    P
      0.8   1.0
O    S
      7.0   0.4
      1.1   0.7
O    P
      1.9   0.5
      0.45  0.6
O    D
      1.0   1.0
O    F
      0.8   1.0
END
)";

BasisSet test_basis_set(bool spherical, const Vec3& shift = {0.0, 0.0, 0.0}) {
  Molecule mol;
  mol.atoms = {{8, {0.1, -0.2, 0.3}}, {1, {1.4, 0.9, -0.5}}};
  for (Atom& atom : mol.atoms)
    for (std::size_t k = 0; k < 3; ++k) atom.r[k] += shift[k];
  std::istringstream in(test_basis);
  return BasisSet::build(mol, BasisLibrary::parse_nwchem(in, "test"), spherical);
}

}  // namespace

// Two normalised s Gaussians squared are spherical Gaussian charge
// distributions of unit charge. Far apart they repel like point charges:
// (aa|bb) -> 1/R. This fixes the overall prefactor of the integral formula
// by electrostatics alone.
TEST_CASE("(ss|ss): two distant charge distributions repel like point charges", "[integrals][eri]") {
  const double R = 30.0;
  const BasisSet basis = BasisSet::from_shells(
      "s", true, {s_shell(0, {0.0, 0.0, 0.0}, 0.7), s_shell(1, {0.0, 0.0, R}, 1.9)});
  const auto& shells = basis.shells();
  EriEngine engine(0);
  const ShellPair aa = make_shell_pair(shells[0], shells[0]);
  const ShellPair bb = make_shell_pair(shells[1], shells[1]);
  CHECK_THAT(engine.compute(aa, bb)[0], WithinRel(1.0 / R, 1e-13));
  CHECK_THAT(engine.compute(bb, aa)[0], WithinRel(1.0 / R, 1e-13));
}

// Closed forms for s functions. Coulomb energy of two spherical Gaussian
// charge distributions with exponents 2a and 2b at distance R:
// erf(sqrt(g) R) / R with g = 2a 2b / (2a + 2b); for R -> 0 this becomes
// 2 sqrt(g / pi).
TEST_CASE("(ss|ss): closed forms", "[integrals][eri]") {
  const double a = 0.9, b = 1.6;
  const double g = 4.0 * a * b / (2.0 * a + 2.0 * b);
  EriEngine engine(0);
  for (double R : {0.0, 0.3, 1.1, 4.0}) {
    INFO("R = " << R);
    const BasisSet basis = BasisSet::from_shells(
        "s", true, {s_shell(0, {0.2, 0.0, 0.0}, a), s_shell(1, {0.2, R, 0.0}, b)});
    const auto& shells = basis.shells();
    const double value =
        engine.compute(make_shell_pair(shells[0], shells[0]), make_shell_pair(shells[1], shells[1]))[0];
    const double exact = R == 0.0 ? 2.0 * std::sqrt(g / M_PI) : std::erf(std::sqrt(g) * R) / R;
    CHECK_THAT(value, WithinRel(exact, 1e-13));
  }
}

TEST_CASE("two-electron tensor: permutational symmetry", "[integrals][eri]") {
  for (const bool spherical : {true, false}) {
    INFO((spherical ? "spherical" : "cartesian"));
    const testing::FullEri eri = testing::full_eri_tensor(test_basis_set(spherical));
    double worst = 0.0;
    for (int i = 0; i < eri.n; ++i)
      for (int j = 0; j < eri.n; ++j)
        for (int k = 0; k < eri.n; ++k)
          for (int l = 0; l < eri.n; ++l) {
            const double v = eri(i, j, k, l);
            worst = std::max({worst, std::abs(v - eri(j, i, k, l)), std::abs(v - eri(i, j, l, k)),
                              std::abs(v - eri(k, l, i, j))});
          }
    CHECK(worst < 1e-12);
  }
}

TEST_CASE("two-electron tensor: positivity and Schwarz inequality", "[integrals][eri]") {
  const testing::FullEri eri = testing::full_eri_tensor(test_basis_set(true));
  for (int i = 0; i < eri.n; ++i)
    for (int j = 0; j < eri.n; ++j) {
      // (ij|ij) is the self-repulsion of the charge distribution ij.
      CHECK(eri(i, j, i, j) >= -1e-14);
      for (int k = 0; k < eri.n; ++k)
        for (int l = 0; l < eri.n; ++l)
          CHECK(eri(i, j, k, l) * eri(i, j, k, l) <= eri(i, j, i, j) * eri(k, l, k, l) * (1.0 + 1e-10) + 1e-20);
    }
}

TEST_CASE("two-electron tensor: translation invariance", "[integrals][eri]") {
  const testing::FullEri eri = testing::full_eri_tensor(test_basis_set(true));
  const testing::FullEri moved = testing::full_eri_tensor(test_basis_set(true, {1.7, -0.6, 2.3}));
  REQUIRE(eri.values.size() == moved.values.size());
  double worst = 0.0;
  for (std::size_t k = 0; k < eri.values.size(); ++k)
    worst = std::max(worst, std::abs(eri.values[k] - moved.values[k]));
  CHECK(worst < 1e-12);
}

TEST_CASE("EriEngine rejects shells beyond its angular momentum", "[integrals][eri]") {
  const BasisSet basis = test_basis_set(true);
  const auto& shells = basis.shells();
  EriEngine engine(1);
  const ShellPair pp = make_shell_pair(shells[1], shells[1]);  // oxygen p
  const ShellPair ff = make_shell_pair(shells[3], shells[3]);  // oxygen f
  CHECK_NOTHROW(engine.compute(pp, pp));
  CHECK_THROWS_AS(engine.compute(ff, pp), std::invalid_argument);
  CHECK_THROWS_AS(engine.compute(pp, ff), std::invalid_argument);
}
