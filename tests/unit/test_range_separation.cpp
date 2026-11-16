// erf-attenuated integrals and range-separated exchange: checks without
// reference code, including the limits named in the specification
// (section 5): omega -> infinity gives the Coulomb integrals, omega -> 0 zero.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "core/xyz.hpp"
#include "df/df_integrals.hpp"
#include "df/df_jk.hpp"
#include "integrals/shell_pair.hpp"
#include "integrals/two_electron.hpp"
#include "scf/rhf.hpp"
#include "xc/xc_evaluator.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

Molecule water() { return read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2o.xyz"); }

Eigen::MatrixXd general_matrix(Eigen::Index n, unsigned seed) {
  Eigen::MatrixXd m(n, n);
  unsigned state = seed;
  for (Eigen::Index i = 0; i < n; ++i)
    for (Eigen::Index j = 0; j < n; ++j) {
      state = state * 1664525u + 1013904223u;
      m(i, j) = static_cast<double>(state >> 8) / static_cast<double>(1u << 24) - 0.5;
    }
  return m;
}

// Largest element of all shell quartet blocks, and largest deviation from the
// Coulomb blocks, for a given omega.
struct Deviation {
  double largest = 0.0, from_coulomb = 0.0;
};
Deviation scan(const BasisSet& basis, double omega) {
  const auto& shells = basis.shells();
  EriEngine attenuated(basis.max_l()), coulomb(basis.max_l());
  Deviation d;
  for (std::size_t a = 0; a < shells.size(); ++a)
    for (std::size_t b = 0; b <= a; ++b) {
      const ShellPair bra = make_shell_pair(shells[a], shells[b]);
      for (std::size_t c = 0; c <= a; ++c) {
        const ShellPair ket = make_shell_pair(shells[c], shells[c]);
        const std::vector<double>& x = attenuated.compute(bra, ket, omega);
        const std::vector<double>& y = coulomb.compute(bra, ket);
        for (std::size_t k = 0; k < x.size(); ++k) {
          d.largest = std::max(d.largest, std::abs(x[k]));
          d.from_coulomb = std::max(d.from_coulomb, std::abs(x[k] - y[k]));
        }
      }
    }
  return d;
}

}  // namespace

TEST_CASE("erf-attenuated integrals: limits of the range parameter", "[integrals][erf][m4]") {
  const Molecule mol = water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("6-31G*"));

  // omega -> infinity: erf(omega r)/r -> 1/r.
  CHECK(scan(basis, 1e8).from_coulomb < 1e-12);
  CHECK(scan(basis, 1e4).from_coulomb < 1e-6);
  // omega -> 0: erf(omega r)/r -> 2 omega / sqrt(pi) -> 0, linearly in omega.
  const double small = scan(basis, 1e-6).largest, smaller = scan(basis, 1e-9).largest;
  CHECK(small < 1e-5);
  CHECK(smaller < 1e-8);
  CHECK_THAT(small / smaller, WithinRel(1e3, 1e-6));
  // In between the integrals lie between 0 and the Coulomb values in size.
  const Deviation mid = scan(basis, 0.33);
  CHECK(mid.largest > 0.1);
  CHECK(mid.from_coulomb > 0.1);

  // Two-centre integrals: same limits.
  const BasisSet aux = BasisSet::build(mol, BasisLibrary::load("def2-universal-JFIT"));
  const Eigen::MatrixXd V = two_center_coulomb(aux);
  CHECK((two_center_coulomb(aux, 1e8) - V).cwiseAbs().maxCoeff() < 1e-10 * V.cwiseAbs().maxCoeff());
  CHECK(two_center_coulomb(aux, 1e-9).cwiseAbs().maxCoeff() < 1e-6);
}

// Closed form for two normalised s Gaussians (exponents a, b) at distance R:
// the charge distributions have exponents 2a, 2b, and
//   (aa|bb)_erf = erf(sqrt(g) R) / R,  1/g = 1/(2a) + 1/(2b) + 1/omega^2.
// (Attenuation adds a third Gaussian width; for omega -> infinity this is the
// Coulomb formula tested in test_two_electron.cpp.)
TEST_CASE("erf-attenuated (ss|ss): closed form", "[integrals][erf]") {
  const double a = 0.9, b = 1.6, omega = 0.4;
  EriEngine engine(0);
  auto s_shell = [](const Vec3& center, double exponent) {
    Shell shell;
    shell.center = center;
    shell.exponents = {exponent};
    shell.coefficients = {1.0};
    return shell;
  };
  for (double R : {0.0, 0.7, 2.5}) {
    INFO("R = " << R);
    const BasisSet basis = BasisSet::from_shells("s", true, {s_shell({0.0, 0.0, 0.0}, a), s_shell({0.0, 0.0, R}, b)});
    const auto& shells = basis.shells();
    const double value = engine.compute(make_shell_pair(shells[0], shells[0]),
                                        make_shell_pair(shells[1], shells[1]), omega)[0];
    const double g = 1.0 / (1.0 / (2.0 * a) + 1.0 / (2.0 * b) + 1.0 / (omega * omega));
    const double exact = R == 0.0 ? 2.0 * std::sqrt(g / M_PI) : std::erf(std::sqrt(g) * R) / R;
    CHECK_THAT(value, WithinRel(exact, 1e-13));
  }
}

TEST_CASE("long-range exchange: limits and both builders", "[scf][jk][erf]") {
  const Molecule mol = water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("6-31G*"));
  const Eigen::MatrixXd D = general_matrix(basis.n_functions(), 9);

  DirectJK direct(basis, DirectJKOptions{0.0});
  std::vector<Eigen::MatrixXd> J, K, K_lr;
  direct.build({D}, &J, &K);
  direct.build_long_range_exchange({D}, 1e8, K_lr);
  CHECK((K_lr[0] - K[0]).cwiseAbs().maxCoeff() < 1e-10);  // omega -> infinity
  direct.build_long_range_exchange({D}, 1e-9, K_lr);
  CHECK(K_lr[0].cwiseAbs().maxCoeff() < 1e-6);            // omega -> 0
  CHECK_THROWS_AS(direct.build_long_range_exchange({D}, 0.0, K_lr), std::invalid_argument);

  // Density fitting: close to the direct result, and equal to the full
  // exchange for omega -> infinity.
  const BasisSet aux = BasisSet::build(mol, BasisLibrary::load("def2-universal-JKFIT"));
  DensityFittedJK fitted(basis, aux);
  std::vector<Eigen::MatrixXd> K_direct, K_fit, K_fit_full, K_fit_inf;
  direct.build_long_range_exchange({D}, 0.33, K_direct);
  fitted.build_long_range_exchange({D}, 0.33, K_fit);
  // The attenuated kernel is smooth, so the fit is very accurate here.
  CHECK((K_fit[0] - K_direct[0]).cwiseAbs().maxCoeff() < 1e-6);
  fitted.build({D}, nullptr, &K_fit_full);
  fitted.build_long_range_exchange({D}, 1e8, K_fit_inf);
  CHECK((K_fit_inf[0] - K_fit_full[0]).cwiseAbs().maxCoeff() < 1e-8);
  CHECK_THROWS_AS(fitted.build_long_range_exchange({D}, -1.0, K_fit), std::invalid_argument);
}

// Specification, section 5: CAM-B3LYP is 0.19 K + 0.46 K_LR with omega = 0.33;
// libxc returns the coefficients in a different split. The conversion:
TEST_CASE("exact-exchange coefficients of the hybrids", "[xc][m4]") {
  const Functional cam = Functional::by_name("CAM-B3LYP");
  CHECK_THAT(cam.exact_exchange, WithinAbs(0.19, 1e-12));
  CHECK_THAT(cam.long_range_exchange, WithinAbs(0.46, 1e-12));
  CHECK_THAT(cam.omega, WithinAbs(0.33, 1e-12));
  CHECK(cam.components.size() == 1);
  CHECK(cam.components[0].libxc_name == "HYB_GGA_XC_CAM_B3LYP");

  const Functional pbe0 = Functional::by_name("pbe0");
  CHECK_THAT(pbe0.exact_exchange, WithinAbs(0.25, 1e-12));
  CHECK(pbe0.long_range_exchange == 0.0);
  CHECK(pbe0.omega == 0.0);

  const Functional b3lyp = Functional::by_name("B3LYP"), b3lyp5 = Functional::by_name("B3LYP5");
  CHECK_THAT(b3lyp.exact_exchange, WithinAbs(0.20, 1e-12));
  CHECK(b3lyp.components[0].libxc_name == "HYB_GGA_XC_B3LYP");    // VWN-RPA
  CHECK(b3lyp5.components[0].libxc_name == "HYB_GGA_XC_B3LYP5");  // VWN5
  // The two B3LYP variants differ in the correlation energy.
  const auto x = make_libxc_evaluator(b3lyp), y = make_libxc_evaluator(b3lyp5);
  const double rho = 0.3, sigma = 0.05;
  double e1, e2, v, vs;
  x->evaluate(1, &rho, &sigma, &e1, &v, &vs);
  y->evaluate(1, &rho, &sigma, &e2, &v, &vs);
  CHECK(std::abs(e1 - e2) > 1e-4);

  // Pure functionals have no exact exchange.
  const HybridParameters pure = libxc_hybrid_parameters("GGA_X_PBE");
  CHECK(pure.exact_exchange == 0.0);
  CHECK(pure.long_range_exchange == 0.0);
  CHECK_THROWS_AS(libxc_hybrid_parameters("NOT_A_FUNCTIONAL"), std::invalid_argument);
}

// The long-range exchange with omega -> infinity is the full exchange, so a
// functional written as c K_LR(omega = 1e8) must give the Hartree-Fock energy.
TEST_CASE("SCF with long-range exchange reduces to Hartree-Fock for large omega", "[scf][erf]") {
  const Molecule mol = water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("STO-3G"));
  Method method;
  method.functional = Functional::by_name("HF");
  method.functional.exact_exchange = 0.3;
  method.functional.long_range_exchange = 0.7;
  method.functional.omega = 1e8;
  const ScfResult split = run_scf(mol, basis, method);
  const ScfResult hf = run_rhf(mol, basis);
  REQUIRE(split.converged);
  CHECK_THAT(split.energy_total, WithinAbs(hf.energy_total, 1e-9));
}
