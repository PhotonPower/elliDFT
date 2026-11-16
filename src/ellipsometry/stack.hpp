#pragma once

// Forward model of ellipsometry: reflection of a plane wave from a stack of
// homogeneous, isotropic layers on a substrate.

#include <memory>
#include <vector>

#include "ellipsometry/materials.hpp"

namespace ellidft::ellipsometry {

// A layer is isotropic, or uniaxial with the optic axis along the surface
// normal: `material` then gives the ordinary index (field in the plane of the
// film), `extraordinary` the index for the field along the normal.
struct Layer {
  std::shared_ptr<Material> material;
  double thickness_nm = 0.0;
  std::shared_ptr<Material> extraordinary;   // empty: isotropic

  Layer() = default;
  Layer(std::shared_ptr<Material> material_, double thickness_nm_, std::shared_ptr<Material> extraordinary_ = nullptr)
      : material(std::move(material_)), thickness_nm(thickness_nm_), extraordinary(std::move(extraordinary_)) {}
};

// Ambient (semi-infinite, where the light comes from), layers from top to
// bottom, substrate.
//
// The substrate is semi-infinite unless substrate_thickness_nm > 0. Then it
// is a thick plate with the medium `backside` (empty: the ambient) behind
// it, and the light reflected from its back face is added incoherently (the
// plate is much thicker than the coherence length): see observe().
struct Stack {
  std::shared_ptr<Material> ambient;
  std::vector<Layer> layers;
  std::shared_ptr<Material> substrate;
  double substrate_thickness_nm = 0.0;
  std::shared_ptr<Material> backside;

  Stack() = default;
  Stack(std::shared_ptr<Material> ambient_, std::vector<Layer> layers_, std::shared_ptr<Material> substrate_)
      : ambient(std::move(ambient_)), layers(std::move(layers_)), substrate(std::move(substrate_)) {}
};

// Amplitude reflection coefficients of the whole stack for s and p
// polarisation. Sign convention of the Fresnel coefficients as in Wikipedia,
// "Fresnel equations" (non-magnetic media): at normal incidence r_p = -r_s.
struct Reflection {
  Complex rs, rp;
};
Reflection reflect(const Stack& stack, double wavelength_nm, double angle_deg);

// Ellipsometric angles: rho = r_p / r_s = tan(Psi) exp(i Delta)
// (Wikipedia, "Ellipsometry"). Psi in [0, 90] degrees, Delta in [0, 360).
// With the sign convention above, a bare transparent substrate has
// Delta = 180 degrees below the Brewster angle and 0 above. Instruments and
// programs that use the opposite sign for r_p report Delta shifted by 180.
struct PsiDelta {
  double psi_deg = 0.0, delta_deg = 0.0;
};
PsiDelta psi_delta(const Reflection& reflection);

// What an ellipsometer sees, including the back face of a thick substrate.
//
// Without a back face: from the coherent reflection, as above.
//
// With a back face (derived here): the beams that leave after m = 1, 2, ...
// round trips in the plate have the amplitudes
//   b_m = t t' r_b a^2 (r' r_b a^2)^(m-1),
// t, t' the transmission of the layer stack downwards and upwards, r' its
// reflection from below, r_b the reflection at the back face, a the
// attenuation of one pass. Beams of different m do not interfere, the s and p
// parts of the same beam do. With
//   I_ss = |r_s|^2 + sum |b_m,s|^2,  I_pp likewise,  I_ps = r_p r_s^* + sum b_m,p b_m,s^*
// (geometric series) the quantities
//   N = (I_ss - I_pp) / (I_ss + I_pp),  C = 2 Re I_ps / (I_ss + I_pp),  S = 2 Im I_ps / (I_ss + I_pp)
// reduce to cos 2Psi, sin 2Psi cos Delta, sin 2Psi sin Delta for a single
// beam. Psi and Delta are DEFINED here as Psi = acos(N) / 2 and
// Delta = atan2(S, C); `polarisation` = sqrt(N^2 + C^2 + S^2) is below one
// when the back face depolarises. How a given instrument reports Psi and
// Delta for partly depolarised light has to be checked against its manual.
struct Observation {
  double N = 0.0, C = 0.0, S = 0.0;
  double psi_deg = 0.0, delta_deg = 0.0;
  double polarisation = 1.0;
  double Rs = 0.0, Rp = 0.0;   // reflected fraction of the intensity
};
Observation observe(const Stack& stack, double wavelength_nm, double angle_deg);

// Psi and Delta of observe().
PsiDelta psi_delta(const Stack& stack, double wavelength_nm, double angle_deg);

}  // namespace ellidft::ellipsometry
