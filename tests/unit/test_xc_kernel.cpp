// Exchange-correlation kernel. Specification, section 7: kernel against the
// numerical derivative of the exchange-correlation matrix.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "core/xyz.hpp"
#include "scf/rhf.hpp"
#include "scf/xc_matrix.hpp"

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

}  // namespace

TEST_CASE("second derivatives of the functionals against finite differences", "[xc][kernel]") {
  for (const char* name : {"LDA", "PBE", "CAM-B3LYP"}) {
    INFO(name);
    const auto xc = make_libxc_evaluator(Functional::by_name(name));
    for (double rho : {0.02, 0.4, 5.0}) {
      for (double s : {0.2, 1.5}) {
        const double sigma = s * std::pow(rho, 8.0 / 3.0);
        INFO("rho = " << rho << ", sigma = " << sigma);
        double v, vs, f_rr, f_rs, f_ss, e, d1, d2, d3;
        xc->evaluate_kernel(1, &rho, &sigma, &v, &vs, &f_rr, &f_rs, &f_ss);
        // The first derivatives agree with evaluate().
        double v_check, vs_check;
        xc->evaluate(1, &rho, &sigma, &e, &v_check, &vs_check);
        CHECK_THAT(v, WithinRel(v_check, 1e-13));
        CHECK_THAT(vs, WithinAbs(vs_check, 1e-13 * std::abs(vs_check) + 1e-300));

        const double h = 1e-4 * rho;
        double vp, vm, vsp, vsm;
        const double rp = rho + h, rm = rho - h;
        xc->evaluate_kernel(1, &rp, &sigma, &vp, &vsp, &d1, &d2, &d3);
        xc->evaluate_kernel(1, &rm, &sigma, &vm, &vsm, &d1, &d2, &d3);
        CHECK_THAT((vp - vm) / (2.0 * h), WithinRel(f_rr, 1e-6));
        if (xc->needs_gradient()) {
          CHECK_THAT((vsp - vsm) / (2.0 * h), WithinAbs(f_rs, 1e-6 * std::abs(f_rs) + 4e-16 * std::abs(vs) / h));
          const double hs = 1e-4 * sigma, sp = sigma + hs, sm = sigma - hs;
          xc->evaluate_kernel(1, &rho, &sp, &vp, &vsp, &d1, &d2, &d3);
          xc->evaluate_kernel(1, &rho, &sm, &vm, &vsm, &d1, &d2, &d3);
          CHECK_THAT((vsp - vsm) / (2.0 * hs), WithinAbs(f_ss, 1e-5 * std::abs(f_ss) + 4e-16 * std::abs(vs) / hs));
          CHECK_THAT((vp - vm) / (2.0 * hs), WithinAbs(f_rs, 1e-5 * std::abs(f_rs) + 4e-16 * std::abs(v) / hs));
        } else {
          CHECK(f_rs == 0.0);
          CHECK(f_ss == 0.0);
        }
      }
    }
  }
  const auto own = make_own_lda_evaluator();
  const double rho = 0.3;
  double a, b, c, d, e;
  CHECK_THROWS_AS(own->evaluate_kernel(1, &rho, nullptr, &a, &b, &c, &d, &e), std::logic_error);
}

TEST_CASE("XC kernel is the derivative of the XC matrix", "[scf][kernel][m5]") {
  const Molecule mol = water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("6-31G*"));
  const MolecularGrid grid = build_molecular_grid(mol, {30, 110});
  const Eigen::MatrixXd D = run_rhf(mol, basis).density;
  const Eigen::Index n = basis.n_functions();
  const Eigen::MatrixXd P = general_matrix(n, 21);                  // unsymmetric trial density
  const Eigen::MatrixXd P_sym = 0.5 * (P + P.transpose());

  for (const char* name : {"LDA", "PBE", "PBE0", "CAM-B3LYP"}) {
    INFO(name);
    const auto xc = make_libxc_evaluator(Functional::by_name(name));
    const std::vector<Eigen::MatrixXd> kernel = integrate_xc_kernel(basis, grid, *xc, D, {P, P_sym});
    REQUIRE(kernel.size() == 2);
    CHECK((kernel[0] - kernel[0].transpose()).cwiseAbs().maxCoeff() < 1e-13);
    // Only the symmetric part of the trial density matters.
    CHECK((kernel[0] - kernel[1]).cwiseAbs().maxCoeff() < 1e-13);

    // Central differences of the XC matrix with Richardson extrapolation.
    auto difference = [&](double h) {
      return Eigen::MatrixXd((integrate_xc(basis, grid, *xc, D + h * P_sym).matrix -
                              integrate_xc(basis, grid, *xc, D - h * P_sym).matrix) / (2.0 * h));
    };
    const double h = 2e-4;
    const Eigen::MatrixXd numeric = (4.0 * difference(h) - difference(2.0 * h)) / 3.0;
    const double scale = kernel[0].cwiseAbs().maxCoeff();
    INFO("largest element " << scale << ", deviation " << (numeric - kernel[0]).cwiseAbs().maxCoeff());
    CHECK((numeric - kernel[0]).cwiseAbs().maxCoeff() < 1e-7 * scale);
  }

  const auto pbe = make_libxc_evaluator(Functional::by_name("PBE"));
  CHECK_THROWS_AS(integrate_xc_kernel(basis, grid, *pbe, D, {Eigen::MatrixXd::Zero(2, 2)}), std::invalid_argument);
  CHECK(integrate_xc_kernel(basis, grid, *pbe, D, {}).empty());
}
