// Linear response properties: checks without reference code
// (specification, sections 6 M5 and 7).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include "core/xyz.hpp"
#include "explicit_ab.hpp"
#include "response/excitations.hpp"
#include "response/polarizability.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;

namespace {

// Water turned into a general orientation, so that no tensor component
// vanishes by symmetry.
Molecule tilted_water() {
  Molecule mol = read_xyz_file(std::string(ELLIDFT_SOURCE_DIR) + "/data/geometries/h2o.xyz");
  const Eigen::Vector3d u = Eigen::Vector3d(0.4, -0.9, 0.6).normalized();
  Eigen::Matrix3d cross;
  cross << 0.0, -u.z(), u.y(), u.z(), 0.0, -u.x(), -u.y(), u.x(), 0.0;
  const Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + std::sin(0.77) * cross + (1.0 - std::cos(0.77)) * cross * cross;
  for (Atom& atom : mol.atoms) {
    const Eigen::Vector3d r = R * Eigen::Vector3d(atom.r[0], atom.r[1], atom.r[2]);
    atom.r = {r.x(), r.y(), r.z()};
  }
  return mol;
}

ScfOptions tight() {
  ScfOptions options;
  options.energy_tolerance = 1e-12;
  options.error_tolerance = 1e-10;
  return options;
}

LinearResponseOptions tight_response() {
  LinearResponseOptions options;
  options.tolerance = 1e-11;
  return options;
}

// All excitations from the explicit matrices:
//   (A-B)^(1/2) (A+B) (A-B)^(1/2) y = omega^2 y,  X+Y = (A-B)^(1/2) y / sqrt(omega),
// transition dipoles t = sqrt(2) mu^T (X+Y).
struct FullSpectrum {
  Eigen::VectorXd omega;
  Eigen::MatrixXd t;  // 3 x n_states
};
FullSpectrum full_spectrum(const testing::ExplicitAB& ab, const Eigen::MatrixXd& mu) {
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> minus(ab.A - ab.B);
  const Eigen::MatrixXd root = minus.eigenvectors() * minus.eigenvalues().cwiseSqrt().asDiagonal() *
                               minus.eigenvectors().transpose();
  const Eigen::MatrixXd M = root * (ab.A + ab.B) * root;
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(0.5 * (M + M.transpose()));
  FullSpectrum s;
  s.omega = eigen.eigenvalues().cwiseSqrt();
  const Eigen::MatrixXd Z = root * eigen.eigenvectors() * s.omega.cwiseSqrt().cwiseInverse().asDiagonal();
  s.t = std::sqrt(2.0) * mu.transpose() * Z;
  return s;
}

}  // namespace

// Acceptance criterion M5: static polarisability from CPKS against finite
// fields (central differences at two field strengths, Richardson
// extrapolation), relative deviation at most 1e-4 for every tensor component.
TEST_CASE("static polarisability: CPKS equals finite fields", "[response][m5]") {
  const Molecule mol = tilted_water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("STO-3G"));
  for (const char* name : {"HF", "PBE", "PBE0", "CAM-B3LYP"}) {
    DYNAMIC_SECTION(name) {
      Method method;
      method.functional = Functional::by_name(name);
      method.grid = {30, 110};
      DirectJK jk(basis, DirectJKOptions{0.0});
      const ScfResult scf = run_scf(mol, basis, method, jk, tight());
      REQUIRE(scf.converged);
      const ResponseOperator op(mol, basis, method, scf, jk);
      const Polarizability cpks = polarizability(op, 0.0, 0.0, tight_response());
      REQUIRE(cpks.converged);
      CHECK(cpks.imag.cwiseAbs().maxCoeff() == 0.0);

      // alpha_kl = d mu_k / d F_l
      Eigen::Matrix3d finite;
      const double h = 2e-3;
      for (int l = 0; l < 3; ++l) {
        auto dipole_at = [&](double strength) {
          ScfOptions options = tight();
          options.electric_field[static_cast<std::size_t>(l)] = strength;
          const ScfResult r = run_scf(mol, basis, method, jk, options);
          REQUIRE(r.converged);
          return Eigen::Vector3d(r.dipole[0], r.dipole[1], r.dipole[2]);
        };
        const Eigen::Vector3d d1 = (dipole_at(h) - dipole_at(-h)) / (2.0 * h);
        const Eigen::Vector3d d2 = (dipole_at(2.0 * h) - dipole_at(-2.0 * h)) / (4.0 * h);
        finite.col(l) = (4.0 * d1 - d2) / 3.0;
      }
      double worst = 0.0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l) {
          INFO("component " << k << l << ": CPKS " << cpks.real(k, l) << ", finite field " << finite(k, l));
          const double relative = std::abs(cpks.real(k, l) - finite(k, l)) / std::abs(finite(k, l));
          worst = std::max(worst, relative);
          CHECK(relative <= 1e-4);
        }
      SUCCEED(name << ": largest relative deviation CPKS vs finite field " << worst);

      // Symmetric and positive definite (specification, section 7).
      CHECK((cpks.real - cpks.real.transpose()).cwiseAbs().maxCoeff() < 1e-8);
      CHECK(Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(cpks.real).eigenvalues().minCoeff() > 0.0);
    }
  }
}

// Acceptance criterion M5: small basis with full diagonalisation; the
// directly solved dynamic polarisability agrees with the sum over all states
// to a relative 1e-8.
TEST_CASE("dynamic polarisability equals the sum over all states", "[response][m5]") {
  const Molecule mol = tilted_water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("STO-3G"));
  for (const char* name : {"HF", "PBE", "PBE0", "CAM-B3LYP"}) {
    DYNAMIC_SECTION(name) {
      Method method;
      method.functional = Functional::by_name(name);
      method.grid = {30, 110};
      DirectJK jk(basis, DirectJKOptions{0.0});
      const ScfResult scf = run_scf(mol, basis, method, jk, tight());
      REQUIRE(scf.converged);
      const ResponseOperator op(mol, basis, method, scf, jk);
      const FullSpectrum spectrum = full_spectrum(testing::explicit_ab(mol, basis, method, scf), op.dipole_vectors());
      REQUIRE(spectrum.omega.minCoeff() > 0.3);  // the frequencies below are off resonance

      double worst = 0.0;
      struct Frequency { double omega, gamma; };
      for (const Frequency& f : {Frequency{0.0, 0.0}, Frequency{0.08, 0.0}, Frequency{0.2, 0.0},
                                 Frequency{0.1, 0.02}, Frequency{0.45, 0.01}}) {
        INFO("omega = " << f.omega << ", gamma = " << f.gamma);
        const Polarizability direct = polarizability(op, f.omega, f.gamma, tight_response());
        REQUIRE(direct.converged);
        // alpha_kl = sum_I 2 omega_I t^k_I t^l_I / (omega_I^2 - z^2), z = omega + i gamma
        const std::complex<double> z(f.omega, f.gamma);
        Eigen::Matrix3cd sum = Eigen::Matrix3cd::Zero();
        for (Eigen::Index I = 0; I < spectrum.omega.size(); ++I)
          sum += (2.0 * spectrum.omega(I) / (spectrum.omega(I) * spectrum.omega(I) - z * z)) *
                 (spectrum.t.col(I) * spectrum.t.col(I).transpose()).cast<std::complex<double>>();
        const double scale = sum.cwiseAbs().maxCoeff();
        const double deviation = std::max((direct.real - sum.real()).cwiseAbs().maxCoeff(),
                                          (direct.imag - sum.imag()).cwiseAbs().maxCoeff()) / scale;
        INFO("relative deviation " << deviation);
        CHECK(deviation <= 1e-8);
        worst = std::max(worst, deviation);
        // Sign: for positive frequency and damping the imaginary part is not negative.
        if (f.gamma > 0.0)
          CHECK(Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(direct.imag).eigenvalues().minCoeff() >= -1e-12);
        else
          CHECK(direct.imag.cwiseAbs().maxCoeff() == 0.0);
      }
      SUCCEED(name << ": largest relative deviation direct vs pole sum " << worst);
    }
  }
}

// Specification, section 7: two distant molecules have the sum of the single
// polarisabilities; atoms are isotropic.
TEST_CASE("polarisability: size consistency and isotropy of atoms", "[response][m5]") {
  const BasisLibrary library = BasisLibrary::load("def2-SVP");
  Method method;
  method.functional = Functional::by_name("PBE0");
  method.grid = {40, 110};

  auto alpha_of = [&](const Molecule& mol) {
    const BasisSet basis = BasisSet::build(mol, library);
    DirectJK jk(basis);
    const ScfResult scf = run_scf(mol, basis, method, jk, tight());
    REQUIRE(scf.converged);
    const ResponseOperator op(mol, basis, method, scf, jk);
    const Polarizability alpha = polarizability(op, 0.0, 0.0, tight_response());
    REQUIRE(alpha.converged);
    return alpha.real;
  };

  SECTION("two hydrogen molecules 150 Bohr apart") {
    Molecule one, two;
    one.atoms = {{1, {0.0, 0.0, 0.0}}, {1, {0.0, 0.0, 1.4}}};
    two.atoms = one.atoms;
    // The second molecule is tilted, so that the sum is not a multiple of one tensor.
    two.atoms.push_back({1, {150.0, 0.0, 0.0}});
    two.atoms.push_back({1, {150.0 + 1.4 * std::sin(0.6), 0.0, 1.4 * std::cos(0.6)}});
    Molecule tilted;
    tilted.atoms = {two.atoms[2], two.atoms[3]};
    const Eigen::Matrix3d a = alpha_of(one), b = alpha_of(tilted), ab = alpha_of(two);
    const double scale = ab.cwiseAbs().maxCoeff();
    INFO("deviation from additivity " << (ab - a - b).cwiseAbs().maxCoeff() << " of " << scale);
    // Dipole-dipole coupling at 150 Bohr is of relative order alpha / R^3 = 1e-6.
    CHECK((ab - a - b).cwiseAbs().maxCoeff() < 1e-5 * scale);
  }
  SECTION("neon is isotropic") {
    Molecule neon;
    neon.atoms = {{10, {0.2, -0.4, 0.1}}};
    const Eigen::Matrix3d alpha = alpha_of(neon);
    const double mean = alpha.trace() / 3.0;
    CHECK(mean > 0.0);
    // The grid is not spherically symmetric, so the tolerance is that of the quadrature.
    CHECK((alpha - mean * Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() < 1e-6 * mean);
  }
}

TEST_CASE("Davidson excitations equal the full diagonalisation", "[response][excitations]") {
  const Molecule mol = tilted_water();
  const BasisSet basis = BasisSet::build(mol, BasisLibrary::load("6-31G*"));
  for (const char* name : {"PBE", "CAM-B3LYP"}) {
    DYNAMIC_SECTION(name) {
      Method method;
      method.functional = Functional::by_name(name);
      method.grid = {20, 50};
      DirectJK jk(basis, DirectJKOptions{0.0});
      const ScfResult scf = run_scf(mol, basis, method, jk, tight());
      REQUIRE(scf.converged);
      const ResponseOperator op(mol, basis, method, scf, jk);
      const testing::ExplicitAB ab = testing::explicit_ab(mol, basis, method, scf);
      const Eigen::MatrixXd mu = op.dipole_vectors();

      ExcitationOptions options;
      options.n_states = 6;
      options.tolerance = 1e-8;

      SECTION("full response") {
        const FullSpectrum exact = full_spectrum(ab, mu);
        const ExcitationResult result = excitations(op, options);
        REQUIRE(result.converged);
        REQUIRE(result.states.size() == 6);
        CHECK_FALSE(result.tamm_dancoff);
        for (std::size_t s = 0; s < 6; ++s) {
          const auto I = static_cast<Eigen::Index>(s);
          INFO("state " << s);
          CHECK_THAT(result.states[s].energy, WithinAbs(exact.omega(I), 1e-9));
          const double f_exact = 2.0 / 3.0 * exact.omega(I) * exact.t.col(I).squaredNorm();
          CHECK_THAT(result.states[s].oscillator_strength, WithinAbs(f_exact, 1e-7));
          // Normalisation (X+Y)^T (X-Y) = 1.
          CHECK_THAT(result.states[s].x_plus_y.dot(result.states[s].x_minus_y), WithinAbs(1.0, 1e-9));
          CHECK(result.states[s].oscillator_strength >= 0.0);
        }
      }
      SECTION("Tamm-Dancoff") {
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> exact(ab.A);
        options.tamm_dancoff = true;
        const ExcitationResult result = excitations(op, options);
        REQUIRE(result.converged);
        CHECK(result.tamm_dancoff);
        for (std::size_t s = 0; s < 6; ++s) {
          const auto I = static_cast<Eigen::Index>(s);
          INFO("state " << s);
          CHECK_THAT(result.states[s].energy, WithinAbs(exact.eigenvalues()(I), 1e-9));
          const Eigen::Vector3d t = std::sqrt(2.0) * mu.transpose() * exact.eigenvectors().col(I);
          CHECK_THAT(result.states[s].oscillator_strength,
                     WithinAbs(2.0 / 3.0 * exact.eigenvalues()(I) * t.squaredNorm(), 1e-7));
          CHECK_THAT(result.states[s].x_plus_y.squaredNorm(), WithinAbs(1.0, 1e-9));
        }
        // Tamm-Dancoff energies lie above the full-response energies.
        options.tamm_dancoff = false;
        const ExcitationResult full = excitations(op, options);
        CHECK(result.states[0].energy >= full.states[0].energy - 1e-10);
      }
    }
  }
}
