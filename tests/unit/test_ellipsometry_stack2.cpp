// Uniaxial layers and the back face of a thick substrate.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <memory>

#include "ellipsometry/analysis.hpp"
#include "ellipsometry/stack.hpp"

using namespace ellidft::ellipsometry;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

const double degree = M_PI / 180.0;
std::shared_ptr<Material> constant(double n, double k = 0.0) { return std::make_shared<ConstantMaterial>(n, k); }

}  // namespace

TEST_CASE("uniaxial layer: limits", "[ellipsometry]") {
  const auto substrate = constant(3.9, 0.02);
  const double wavelength = 600.0, d = 150.0;
  auto stack_with = [&](double n_o, double n_e) {
    return Stack(constant(1.0), {Layer(constant(n_o), d, n_e > 0.0 ? constant(n_e) : nullptr)}, substrate);
  };
  SECTION("equal indices: the isotropic layer") {
    for (double angle : {0.0, 40.0, 70.0}) {
      const Reflection a = reflect(stack_with(1.6, 1.6), wavelength, angle), b = reflect(stack_with(1.6, -1.0), wavelength, angle);
      CHECK_THAT(std::abs(a.rs - b.rs), WithinAbs(0.0, 1e-15));
      CHECK_THAT(std::abs(a.rp - b.rp), WithinAbs(0.0, 1e-15));
    }
  }
  SECTION("normal incidence: the field has no component along the axis") {
    const Reflection a = reflect(stack_with(1.6, 1.9), wavelength, 0.0), b = reflect(stack_with(1.6, -1.0), wavelength, 0.0);
    CHECK_THAT(std::abs(a.rs - b.rs), WithinAbs(0.0, 1e-15));
    CHECK_THAT(std::abs(a.rp - b.rp), WithinAbs(0.0, 1e-15));
  }
  SECTION("s polarisation never sees the extraordinary index, p polarisation does") {
    const Reflection a = reflect(stack_with(1.6, 1.9), wavelength, 65.0), b = reflect(stack_with(1.6, -1.0), wavelength, 65.0);
    CHECK_THAT(std::abs(a.rs - b.rs), WithinAbs(0.0, 1e-15));
    CHECK(std::abs(a.rp - b.rp) > 1e-3);
  }
  SECTION("p polarisation: a layer of phase thickness pi is absent") {
    // normal component of the extraordinary wave: q = (n_o / n_e) sqrt(n_e^2 - sin^2(theta))
    const double n_o = 1.6, n_e = 1.75, angle = 60.0, s = std::sin(angle * degree);
    const double q = n_o / n_e * std::sqrt(n_e * n_e - s * s);
    const Stack stack(constant(1.0), {Layer(constant(n_o), wavelength / (2.0 * q), constant(n_e))}, substrate);
    const Reflection with = reflect(stack, wavelength, angle), bare = reflect(Stack(constant(1.0), {}, substrate), wavelength, angle);
    CHECK_THAT(std::abs(with.rp - bare.rp), WithinAbs(0.0, 1e-14));
    CHECK(std::abs(with.rs - bare.rs) > 1e-3);   // for s the phase thickness is another one
  }
  SECTION("a thick, weakly absorbing uniaxial layer has its Brewster angle at tan^2 = n_e^2 (n_o^2 - 1) / (n_e^2 - 1)") {
    // (from eta_p of the two media being equal; for n_e = n_o the usual tan = n)
    const double n_o = 1.5, n_e = 1.8;
    const double brewster = std::atan(std::sqrt(n_e * n_e * (n_o * n_o - 1.0) / (n_e * n_e - 1.0))) / degree;
    const Stack stack(constant(1.0), {Layer(constant(n_o, 1e-3), 2e6, constant(n_e, 1e-3))}, substrate);
    CHECK(std::abs(reflect(stack, wavelength, brewster).rp) < 2e-3);
    CHECK(std::abs(reflect(stack, wavelength, brewster - 5.0).rp) > 2e-2);
    CHECK(std::abs(reflect(stack, wavelength, std::atan(n_o) / degree).rp) > 1e-2);   // not at the isotropic angle
  }
}

TEST_CASE("back face of a thick substrate", "[ellipsometry]") {
  const double wavelength = 550.0, angle = 60.0;
  SECTION("glass plate in air: R_total = 2 R / (1 + R) for each polarisation") {
    Stack plate(constant(1.0), {}, constant(1.5));
    const Reflection front = reflect(plate, wavelength, angle);
    plate.substrate_thickness_nm = 1e6;
    const Observation o = observe(plate, wavelength, angle);
    const double Rs = std::norm(front.rs), Rp = std::norm(front.rp);
    CHECK_THAT(o.Rs, WithinRel(2.0 * Rs / (1.0 + Rs), 1e-13));
    CHECK_THAT(o.Rp, WithinRel(2.0 * Rp / (1.0 + Rp), 1e-13));
    CHECK(o.polarisation < 1.0);         // the back face depolarises
    CHECK(o.polarisation > 0.9);
    // Psi and Delta as defined from N, C, S
    CHECK_THAT(std::cos(2.0 * o.psi_deg * degree), WithinAbs(o.N, 1e-13));
    const PsiDelta pd = psi_delta(plate, wavelength, angle);
    CHECK(pd.psi_deg == o.psi_deg);
    CHECK(pd.delta_deg == o.delta_deg);
  }
  SECTION("without a back face observe() is the coherent reflection") {
    const Stack stack(constant(1.0), {Layer(constant(1.6), 120.0)}, constant(3.9, 0.02));
    const Observation o = observe(stack, wavelength, angle);
    const PsiDelta pd = psi_delta(reflect(stack, wavelength, angle));
    CHECK_THAT(o.psi_deg, WithinAbs(pd.psi_deg, 1e-10));
    CHECK_THAT(o.delta_deg, WithinAbs(pd.delta_deg, 1e-10));
    CHECK_THAT(o.polarisation, WithinRel(1.0, 1e-13));
  }
  SECTION("an opaque plate and an index-matched back side leave the front reflection") {
    const std::vector<Layer> layers = {Layer(constant(1.6), 120.0)};
    Stack opaque(constant(1.0), layers, constant(1.5, 0.01));
    const PsiDelta front = psi_delta(opaque, wavelength, angle);
    opaque.substrate_thickness_nm = 1e6;
    CHECK_THAT(psi_delta(opaque, wavelength, angle).psi_deg, WithinAbs(front.psi_deg, 1e-9));
    CHECK_THAT(psi_delta(opaque, wavelength, angle).delta_deg, WithinAbs(front.delta_deg, 1e-9));

    Stack matched(constant(1.0), layers, constant(1.5));
    const PsiDelta front_matched = psi_delta(matched, wavelength, angle);
    matched.substrate_thickness_nm = 1e6;
    matched.backside = constant(1.5);
    CHECK_THAT(psi_delta(matched, wavelength, angle).psi_deg, WithinAbs(front_matched.psi_deg, 1e-9));
    CHECK_THAT(observe(matched, wavelength, angle).polarisation, WithinRel(1.0, 1e-12));
    matched.backside = nullptr;   // air behind: the back face changes the result
    CHECK(std::abs(psi_delta(matched, wavelength, angle).psi_deg - front_matched.psi_deg) > 0.1);
  }
  SECTION("a lossless film on a lossless plate reflects at most everything") {
    Stack stack(constant(1.0), {Layer(constant(2.1), 90.0), Layer(constant(1.38), 140.0)}, constant(1.52));
    stack.substrate_thickness_nm = 5e5;
    for (double w = 400.0; w <= 800.0; w += 50.0) {
      const Observation o = observe(stack, w, 55.0);
      CHECK(o.Rs < 1.0);
      CHECK(o.Rp < 1.0);
      CHECK(o.Rs > std::norm(reflect(stack, w, 55.0).rs));   // more than the front alone
    }
  }
}

TEST_CASE("fit of a film on glass needs the back face", "[ellipsometry]") {
  auto make = [](double thickness, bool back_face) {
    auto film = std::make_shared<CauchyMaterial>(1.58, 0.008, 0.0);
    auto stack = std::make_shared<Stack>(constant(1.0), std::vector<Layer>{Layer(film, thickness)}, constant(1.52));
    if (back_face) stack->substrate_thickness_nm = 1e6;
    return stack;
  };
  const auto truth = make(180.0, true);
  std::vector<Measurement> data;
  for (double angle : {55.0, 65.0})
    for (double wavelength = 400.0; wavelength <= 900.0; wavelength += 25.0) {
      const PsiDelta pd = psi_delta(*truth, wavelength, angle);
      Measurement m;
      m.wavelength_nm = wavelength;
      m.angle_deg = angle;
      m.psi_deg = pd.psi_deg;
      m.delta_deg = pd.delta_deg;
      data.push_back(m);
    }
  const auto with = make(150.0, true), without = make(150.0, false);
  const StackFit good = fit_stack(*with, {thickness_parameter(*with, 0, "thickness")}, data, Residuals::PsiDelta);
  const StackFit poor = fit_stack(*without, {thickness_parameter(*without, 0, "thickness")}, data, Residuals::PsiDelta);
  CHECK_THAT(good.fit.parameters[0], WithinAbs(180.0, 1e-5));
  CHECK(good.deviation.rms_psi_deg < 1e-6);
  CHECK(poor.deviation.rms_psi_deg > 0.5);   // the model without back face cannot describe the data
}
