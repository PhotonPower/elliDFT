// Module optics: from polarisabilities to n and k (specification, section 8).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <complex>

#include "optics/optics.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Complex = std::complex<double>;

TEST_CASE("unit conversions of the specification", "[optics]") {
  // N in 1/A^3 = 0.602214076 * rho / M
  CHECK_THAT(optics::number_density(1.05, 104.15), WithinRel(0.602214076 * 1.05 / 104.15, 1e-15));
  // Photon energy in eV = 1239.841984 / wavelength in nm; 1 Hartree = 27.211386245988 eV
  CHECK_THAT(optics::photon_energy_hartree(589.0), WithinRel(1239.841984 / 589.0 / 27.211386245988, 1e-15));
  CHECK_THAT(optics::wavelength_nm(optics::photon_energy_hartree(432.1)), WithinRel(432.1, 1e-14));
  CHECK(optics::angstrom3_per_au == 0.148184711);
  CHECK_THROWS_AS(optics::number_density(0.0, 100.0), std::invalid_argument);
  CHECK_THROWS_AS(optics::photon_energy_hartree(-1.0), std::invalid_argument);
}

// The guide values of the specification (section 8, "Richtwerte für M6"):
// polarisability volumes back-calculated from refractive index and density.
// Polystyrene: M = 104.15, rho = 1.05, n = 1.59 -> 13.3 A^3;
// PMMA: M = 100.12, rho = 1.18, n = 1.49 -> 9.7 A^3 (given to three digits).
TEST_CASE("Lorentz-Lorenz reproduces the guide values of the specification", "[optics]") {
  const double N_ps = optics::number_density(1.05, 104.15);
  const double N_pmma = optics::number_density(1.18, 100.12);
  CHECK_THAT(optics::lorentz_lorenz_inverse(1.59, N_ps).real(), WithinAbs(13.3, 0.05));
  CHECK_THAT(optics::lorentz_lorenz_inverse(1.49, N_pmma).real(), WithinAbs(9.7, 0.05));
  // and back
  CHECK_THAT(optics::lorentz_lorenz(optics::lorentz_lorenz_inverse(1.59, N_ps), N_ps).real(), WithinRel(1.59, 1e-13));
  CHECK_THAT(optics::lorentz_lorenz(13.3, N_ps).real(), WithinAbs(1.59, 0.005));
}

TEST_CASE("Lorentz-Lorenz: limits, absorption and round trip", "[optics]") {
  const double N = 0.006;
  // No polarisability: vacuum.
  CHECK(optics::lorentz_lorenz(0.0, N) == Complex(1.0, 0.0));
  // Dilute limit: n ~ 1 + 2 pi N alpha'.
  CHECK_THAT(optics::lorentz_lorenz(1e-4, N).real(), WithinRel(1.0 + 2.0 * M_PI * N * 1e-4, 1e-9));
  // Absorption: Im alpha' > 0 gives k > 0, and n stays positive.
  const Complex n = optics::lorentz_lorenz(Complex(12.0, 0.8), N);
  CHECK(n.real() > 1.0);
  CHECK(n.imag() > 0.0);
  // The defining equation holds.
  const Complex lhs = (n * n - 1.0) / (n * n + 2.0);
  const Complex rhs = 4.0 * M_PI / 3.0 * N * Complex(12.0, 0.8);
  CHECK_THAT(lhs.real(), WithinRel(rhs.real(), 1e-13));
  CHECK_THAT(lhs.imag(), WithinRel(rhs.imag(), 1e-13));
  // Round trip for complex values.
  const Complex alpha = optics::lorentz_lorenz_inverse(Complex(1.7, 0.05), N);
  const Complex back = optics::lorentz_lorenz(alpha, N);
  CHECK_THAT(back.real(), WithinRel(1.7, 1e-13));
  CHECK_THAT(back.imag(), WithinRel(0.05, 1e-12));
  // Polarisation catastrophe.
  CHECK_THROWS_AS(optics::lorentz_lorenz(3.0 / (4.0 * M_PI * N), N), std::domain_error);
}

TEST_CASE("chain frame and rotation of tensors", "[optics]") {
  // A zig-zag chain along the direction (1, 1, 0)/sqrt(2).
  Molecule chain;
  const Eigen::Vector3d axis = Eigen::Vector3d(1.0, 1.0, 0.0).normalized();
  const Eigen::Vector3d side(0.0, 0.0, 1.0);
  for (int k = 0; k < 9; ++k) {
    const Eigen::Vector3d r = 2.4 * k * axis + (k % 2 == 0 ? 0.8 : -0.8) * side + Eigen::Vector3d(3.0, -1.0, 0.5);
    chain.atoms.push_back({6, {r.x(), r.y(), r.z()}});
  }
  const Eigen::Matrix3d frame = optics::chain_frame(chain);
  CHECK((frame.transpose() * frame - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff() < 1e-13);
  CHECK_THAT(frame.determinant(), WithinAbs(1.0, 1e-13));
  // z of the chain frame is the chain direction (up to sign).
  CHECK_THAT(std::abs(frame.col(2).dot(axis)), WithinAbs(1.0, 1e-12));

  // A tensor that is diagonal along the chain: 10 along, 4 and 3 across.
  const Eigen::Vector3d across1 = side, across2 = axis.cross(side);
  const Eigen::Matrix3d lab = 10.0 * axis * axis.transpose() + 4.0 * across1 * across1.transpose() +
                              3.0 * across2 * across2.transpose();
  const Eigen::Matrix3cd rotated = optics::rotate_into_frame(lab.cast<Complex>(), frame);
  CHECK_THAT(rotated(2, 2).real(), WithinAbs(10.0, 1e-12));
  CHECK_THAT((rotated(0, 0) + rotated(1, 1)).real(), WithinAbs(7.0, 1e-12));
  CHECK_THAT(rotated.trace().real(), WithinAbs(17.0, 1e-12));  // the trace is invariant
  const optics::UniaxialPolarizability u = optics::uniaxial(rotated);
  CHECK_THAT(u.parallel.real(), WithinAbs(10.0, 1e-12));
  CHECK_THAT(u.perpendicular.real(), WithinAbs(3.5, 1e-12));
  CHECK_THAT(u.mean().real(), WithinAbs(17.0 / 3.0, 1e-12));
  CHECK_THAT(u.anisotropy().real(), WithinAbs(6.5, 1e-12));

  Molecule single;
  single.atoms = {{6, {0.0, 0.0, 0.0}}};
  CHECK_THROWS_AS(optics::chain_frame(single), std::invalid_argument);
}

TEST_CASE("increments per repeat unit", "[optics]") {
  // alpha(m) = end groups + m * unit: the increments equal the unit exactly.
  Eigen::Matrix3cd ends = Eigen::Matrix3cd::Zero(), unit = Eigen::Matrix3cd::Zero();
  ends.diagonal() << Complex(2.0), Complex(1.0), Complex(3.0);
  unit.diagonal() << Complex(8.0, 0.1), Complex(7.0, 0.2), Complex(12.0, 0.3);
  std::vector<Eigen::Matrix3cd> oligomers;
  for (int m = 1; m <= 4; ++m) oligomers.push_back(ends + static_cast<double>(m) * unit);
  const auto inc = optics::increments(oligomers);
  REQUIRE(inc.size() == 3);
  for (const auto& tensor : inc) CHECK((tensor - unit).cwiseAbs().maxCoeff() < 1e-13);
  CHECK(optics::increments({}).empty());
  CHECK(optics::increments({unit}).empty());
}

TEST_CASE("orientation average of a uniaxial film", "[optics]") {
  const optics::UniaxialPolarizability u{Complex(15.0, 0.4), Complex(9.0, 0.1)};
  // S = 0: isotropic, both components equal the mean.
  const auto iso = optics::oriented(u, 0.0);
  CHECK(std::abs(iso.ordinary - u.mean()) < 1e-14);
  CHECK(std::abs(iso.extraordinary - u.mean()) < 1e-14);
  // S = 1: chains along the normal; extraordinary = parallel, ordinary = perpendicular.
  const auto aligned = optics::oriented(u, 1.0);
  CHECK(std::abs(aligned.extraordinary - u.parallel) < 1e-14);
  CHECK(std::abs(aligned.ordinary - u.perpendicular) < 1e-14);
  // S = -1/2: chains in the plane; extraordinary = perpendicular.
  const auto in_plane = optics::oriented(u, -0.5);
  CHECK(std::abs(in_plane.extraordinary - u.perpendicular) < 1e-14);
  CHECK(std::abs(in_plane.ordinary - 0.5 * (u.parallel + u.perpendicular)) < 1e-14);
  // The mean (alpha_e + 2 alpha_o)/3 does not depend on S.
  for (double S : {-0.3, 0.2, 0.7}) {
    const auto o = optics::oriented(u, S);
    CHECK(std::abs((o.extraordinary + 2.0 * o.ordinary) / 3.0 - u.mean()) < 1e-14);
  }
  CHECK_THROWS_AS(optics::oriented(u, 1.2), std::invalid_argument);
}

TEST_CASE("optical constants from the polarisability per repeat unit", "[optics]") {
  // A repeat unit with the polystyrene guide value 13.3 A^3, isotropic.
  const double N = optics::number_density(1.05, 104.15);
  const double alpha_au = 13.3 / optics::angstrom3_per_au;
  Eigen::Matrix3cd alpha = Eigen::Matrix3cd::Zero();
  alpha.diagonal().setConstant(Complex(alpha_au, 0.0));
  const optics::OpticalConstants c = optics::optical_constants(589.0, alpha, N);
  CHECK(c.wavelength_nm == 589.0);
  CHECK_THAT(c.isotropic.real(), WithinAbs(1.59, 0.005));
  CHECK(c.isotropic.imag() == 0.0);
  CHECK(std::abs(c.ordinary - c.isotropic) < 1e-14);
  CHECK(std::abs(c.extraordinary - c.isotropic) < 1e-14);

  // Anisotropic unit with absorption: for S > 0 the extraordinary index is
  // the larger one if the chain direction is more polarisable.
  alpha(2, 2) = Complex(1.3 * alpha_au, 0.02 * alpha_au);
  const optics::OpticalConstants film = optics::optical_constants(589.0, alpha, N, 0.4);
  CHECK(film.extraordinary.real() > film.ordinary.real());
  CHECK(film.extraordinary.imag() > 0.0);
  CHECK(film.isotropic.real() > film.ordinary.real());
  CHECK(film.isotropic.real() < film.extraordinary.real());
}
