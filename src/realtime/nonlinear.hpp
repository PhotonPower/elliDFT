#pragma once

// Frequency-dependent non-linear optical response of a molecule from
// real-time propagations in a field that oscillates with one frequency.
//
// Conventions (defined here; the literature uses several that differ by
// numerical factors):
//
//   mu_i = mu0_i + alpha_ij F_j + 1/2 beta_ijk F_j F_k + 1/6 gamma_ijkl F_j F_k F_l
//
// for static fields ("Taylor convention"). The frequency-dependent
// quantities are defined such that all of them become these derivatives for
// omega -> 0. With the field F(t) = F0 + F1 cos(omega t), all vectors along
// fixed directions, collecting powers and harmonics (cos^2 = (1 + cos 2wt)/2,
// cos^3 = (3 cos wt + cos 3wt)/4):
//
//   second order, cos(2wt):  1/4 beta(-2w; w, w) F1 F1            second-harmonic generation
//   second order, cos(wt):       beta(-w; w, 0)  F1 F0            electro-optic (Pockels) effect
//   second order, constant:  1/4 beta(0; w, -w)  F1 F1 + 1/2 beta(0; 0, 0) F0 F0
//   third order,  cos(3wt):  1/24 gamma(-3w; w, w, w) F1^3        third-harmonic generation
//   third order,  cos(wt):   1/8 gamma(-w; w, w, -w) F1^3         optical Kerr effect (intensity-dependent
//                                                                 polarisability); its part in quadrature
//                                                                 (sin wt) is two-photon absorption
//                            + 1/2 gamma(-w; w, 0, 0) F1 F0^2     electro-optic (dc) Kerr effect

#include <array>
#include <complex>
#include <vector>

#include "realtime/propagator.hpp"

namespace ellidft {

// Smooth switch from 0 to 1 on [0, 1]; the first three derivatives vanish at
// both ends (polynomial of degree 7).
double smooth_switch(double x);

struct HarmonicFieldOptions {
  double omega = 0.0;                    // Hartree, > 0
  Vec3 optical{0.0, 0.0, 0.0};           // amplitude of the cos(omega t) part for lambda = 1
  Vec3 static_field{0.0, 0.0, 0.0};      // constant part for lambda = 1
  double amplitude = 0.01;               // a: four propagations with lambda = +-a, +-2a (see dipole_orders)
  double ramp_time = 150.0;              // the field is switched on smoothly during this time
  int n_periods = 2;                     // length of the analysis window after the ramp, in optical periods
  RealTimeOptions propagation;           // time step, propagator, tolerances (n_steps and field are set here)
};

// Harmonic content of the dipole moment, order by order in the field, in the
// window after the switch-on:
//   mu_n(t) = sum_m  cosine[n-1][m] cos(m w t) + sine[n-1][m] sin(m w t),   n = 1, 2, 3;  m = 0 ... 3,
// from a linear least-squares fit. `residual` is the root mean square of what
// the fit leaves, relative to that of the signal: it shows remaining
// transients (switch-on too fast, or a resonance).
struct HarmonicResponse {
  std::array<std::array<Vec3, 4>, 3> cosine{}, sine{};
  std::array<double, 3> residual{};
  std::array<double, 3> signal_rms{};   // root mean square of mu_n in the window (all components)
  int fock_builds = 0;
  bool steps_converged = true;
  double total_time = 0.0;
};
HarmonicResponse harmonic_response(const RealTimePropagator& propagator, const HarmonicFieldOptions& options);

// Two-photon absorption from the energy a pulse leaves in the molecule.
// Pulse: F(t) = F0 sin^2(pi t / T) cos(omega (t - T/2)) e, 0 <= t <= T. The
// absorbed energy is expanded as W(F0) = W2 F0^2 + W4 F0^4 (two propagations,
// F0 = a and 2a); W2 is one-photon absorption (from the spectral wings of the
// pulse, if omega is below the first allowed excitation), W4 two-photon
// absorption.
//
// Cross-section, DEFINED HERE through the number of absorbed photons,
//   dN_photons / dt = sigma2 * (I / omega)^2,
// with the cycle-averaged intensity I(t) = c F(t)^2 / (8 pi) in atomic units:
//   sigma2 = W4 * omega * (8 pi / c)^2 / integral_0^T sin^8(pi t / T) dt.
// It is an average over the two-photon spectrum of the pulse (width about
// 4 pi / T), for molecules of fixed orientation and light polarised along e.
// Other definitions (per excitation instead of per photon, orientational
// average) differ by numerical factors.
struct TwoPhotonOptions {
  double omega = 0.0;              // Hartree
  Vec3 polarization{0.0, 0.0, 1.0};
  double amplitude = 0.01;         // a
  double duration = 300.0;         // T
  // Constant field in addition to the pulse (not scaled with the pulse). For
  // a propagator that was built from the SCF in this same field.
  Vec3 static_field{0.0, 0.0, 0.0};
  RealTimeOptions propagation;     // n_steps and field are set here; Composition4 is advisable (energy)
};
struct TwoPhotonResult {
  double one_photon_coefficient = 0.0;    // W2, Hartree per (field a.u.)^2
  double two_photon_coefficient = 0.0;    // W4, Hartree per (field a.u.)^4
  double cross_section_au = 0.0;          // sigma2 in Bohr^4 * (atomic unit of time)
  double cross_section_gm = 0.0;          // in 1e-50 cm^4 s / photon (Goeppert-Mayer)
  int fock_builds = 0;
  bool steps_converged = true;
};
TwoPhotonResult two_photon_absorption(const RealTimePropagator& propagator, const TwoPhotonOptions& options);

// ---------------------------------------------------------------------------
// Quantities for disordered materials (decisions of the client, 2026-11-15)
// ---------------------------------------------------------------------------

// Seven directions with weights that average any even polynomial of degree 4
// in the direction exactly over the sphere: the three axes (weight 2/15 each)
// and the four space diagonals (3/20 each); these are the 14 points of the
// Lebedev grid of degree 5 with opposite points merged. Both u.gamma.uuu and
// |u.S.u|^2 are such polynomials, so their orientational averages need only
// these directions.
struct AveragingDirection {
  Vec3 u;
  double weight;
};
std::array<AveragingDirection, 7> averaging_directions();

// Orientational average of gamma(-w; w, w, -w) for parallel fields,
//   <u.gamma.uuu> = (1/15) sum_ij (gamma_iijj + gamma_ijij + gamma_ijji),
// from seven harmonic analyses. options.optical and static_field are ignored.
struct KerrAverage {
  double gamma = 0.0;                    // atomic units
  double gamma_quadrature = 0.0;         // part in quadrature (two-photon absorption); zero off resonance
  std::array<double, 7> by_direction{};
  // what the third-order fits leave, relative to the strongest third-order
  // signal among the seven directions (directions without response count little)
  double largest_residual = 0.0;
  int fock_builds = 0;
  bool steps_converged = true;
};
KerrAverage kerr_average(const RealTimePropagator& propagator, const HarmonicFieldOptions& options);

// Intensity-dependent refractive index n = n0 + n2 I of a disordered material
// of molecules with the averaged gamma above, number density N and linear
// index n0, with the Lorentz local field L = (n0^2 + 2) / 3 for every field
// (the approximation the specification prescribes for the linear index):
//
//   delta n = (pi / 4) N L^4 gamma |E|^2 / n0,   I = n0 c |E|^2 / (8 pi)
//   =>  n2 = 2 pi^2 N L^4 gamma / (n0^2 c)        (atomic units).
//
// Derived here from the Lorentz-Lorenz relation with the polarisability
// alpha + gamma F_loc^2 / 8 at the optical frequency (see the table at the top
// of this file); it agrees with n = n0 + 3 chi3 |E|^2 / (8 n0) (Wikipedia,
// "Kerr effect") for chi3 = 4 pi N L^4 gamma / 6. Electronic part only: no
// thermal contribution, no reorientation of molecules.
struct NonlinearIndex {
  double local_field_factor = 0.0;
  double n2_au = 0.0;            // per atomic unit of intensity
  double n2_m2_per_W = 0.0;
  double n2_cm2_per_W = 0.0;
};
NonlinearIndex nonlinear_index(double gamma_au, double number_density_per_A3, double n0);

// The pulse of two_photon_absorption for unit amplitude and its transforms:
//   pulse_transform:             integral_0^T sin^2(pi t/T) cos(w0 (t - T/2)) e^{i w t} dt
//   envelope_squared_transform:  integral_0^T sin^4(pi t/T) e^{i d t} dt
std::complex<double> pulse_transform(double omega0, double duration, double omega);
std::complex<double> envelope_squared_transform(double duration, double delta);

// Two-photon transition strengths of given excited states.
//
// Second-order perturbation theory for the pulse F0 g(t) cos(w t) e with
// 2 w close to the excitation energy w_f gives the amplitude
//   |c_f| = (F0^2 / 8) |e.S_f.e| |integral g^2 e^{i (w_f - 2w) t} dt|,
//   S_f = sum_n [ mu_fn mu_n0 / (w_n - w) + mu_fn mu_n0 / (w_n - w) ]   (two equal photons),
// and the absorbed energy w_f |c_f|^2. For every state one pulse at w_f / 2
// is propagated (two amplitudes, see two_photon_absorption); the F0^4 parts
//   W4_k = sum_f  w_f |e.S_f.e|^2 |G(w_f - 2 w_k)|^2 / 64
// are solved for the strengths |e.S_f.e|^2, which removes the overlap of
// neighbouring states within the width of the pulse. Derived here.
//
// With `average` the strength is the orientational average over the seven
// directions, otherwise that along `polarization`.
struct TwoPhotonStrengthOptions {
  std::vector<double> state_energies;   // Hartree; from the module response
  bool average = true;
  Vec3 polarization{0.0, 0.0, 1.0};
  double amplitude = 0.01;
  double duration = 300.0;
  RealTimeOptions propagation;
};
struct TwoPhotonState {
  double energy = 0.0;      // Hartree
  double strength = 0.0;    // <|e.S.e|^2> or |e.S.e|^2, atomic units
};
struct TwoPhotonStrengths {
  std::vector<TwoPhotonState> states;
  int fock_builds = 0;
  bool steps_converged = true;
};
TwoPhotonStrengths two_photon_strengths(const RealTimePropagator& propagator, const TwoPhotonStrengthOptions& options);

// Cross-section in GM for photons of energy omega, with a Lorentzian line of
// half width `half_width` (Hartree) for every state:
//
//   sigma2(w) = 4 pi^3 alpha^2 w^2 sum_f strength_f g(2w - w_f),
//   g(x) = (1/pi) half_width / (x^2 + half_width^2),
//
// in Bohr^4 times the atomic unit of time; alpha = 1/c. Counted per absorbed
// photon, one beam. Derived here from the rate 2 pi |c_f|^2-per-time and the
// photon flux I / w with I = c F0^2 / (8 pi). It has the form
// N pi^3 alpha a0^5 w^2 / c with N = 4 that a web search reports for
// Beerepoot et al., PCCP 17, 19306 (2015); that paper could not be read.
double two_photon_cross_section_gm(const std::vector<TwoPhotonState>& states, double omega, double half_width);

}  // namespace ellidft
