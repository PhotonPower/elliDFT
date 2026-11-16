#pragma once

// Real-time propagation of the Kohn-Sham (or Hartree-Fock) orbitals of a
// closed-shell molecule after a short electric field pulse ("delta kick").
// The Fourier-Laplace transform of the induced dipole moment gives the
// polarisability alpha(omega + i gamma); see spectrum.hpp.
//
// Purpose in elliDFT: an independent second route to the quantities of the
// module response (the same alpha must come out for weak pulses), and
// broadband spectra from one propagation.

#include <complex>
#include <functional>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "core/molecule.hpp"
#include "grid/molecular_grid.hpp"
#include "scf/jk.hpp"
#include "scf/rhf.hpp"
#include "xc/xc_evaluator.hpp"

namespace ellidft {

enum class Propagator {
  // U = exp(-i dt F_mid) with F_mid = (F_n + F_n+1) / 2. Error of second
  // order in the time step.
  Midpoint,
  // Magnus expansion to fourth order with the Fock matrix at the two
  // Gauss-Legendre points of the step, taken from a polynomial through the
  // Fock matrices of this and earlier steps. Error of fourth order; the
  // cheapest choice for weak pulses, but the energy is not conserved for
  // strong fields and coarse steps.
  Magnus4,
  // Three midpoint steps of lengths w1 dt, w0 dt, w1 dt with
  // w1 = 1 / (2 - 2^(1/3)), w0 = 1 - 2 w1 < 0 ("triple jump"). Error of
  // fourth order; inherits the energy conservation of the midpoint step.
  // About four times the cost of Midpoint or Magnus4 per step (measured).
  Composition4,
};

struct RealTimeOptions {
  double time_step = 0.1;   // atomic units of time (1 a.u. = 24.19 as)
  int n_steps = 1000;
  Propagator propagator = Propagator::Magnus4;

  // The Fock matrix inside a step depends on the orbitals at its end. The
  // step is iterated until no element of the Fock matrix at the end
  // (orthonormal basis, Hartree) changes by more than fock_tolerance; at most
  // max_fock_builds_per_step iterations, one Fock matrix each.
  double fock_tolerance = 1e-9;
  int max_fock_builds_per_step = 8;

  // Homogeneous electric field during the propagation, as a function of the
  // time in atomic units (empty: none). Same convention as
  // ScfOptions::electric_field: H(t) = H0 - mu.F(t), the electrons feel
  // +F(t).r. The function is also called for times slightly outside
  // [0, n_steps * time_step] (Composition4).
  std::function<Vec3(double)> field;
};

struct RealTimeTrajectory {
  double time_step = 0.0;
  // Values at t = 0 (directly after the pulse), dt, 2 dt, ...
  std::vector<Vec3> dipole;      // e * Bohr, nuclei minus electrons, origin (0,0,0)
  std::vector<double> energy;    // Hartree; the energy functional without the interaction with the field

  // Diagnostics.
  int fock_builds = 0;
  int max_fock_builds_in_a_step = 0;
  bool steps_converged = true;             // the iteration of every step reached fock_tolerance
  double orthonormality_error = 0.0;       // largest element of C^+ C - 1 at the end
  double energy_drift = 0.0;               // largest |E(t) - E(0)|; an error measure only without a field
};

class RealTimePropagator {
 public:
  // `scf` must be the converged ground state for `molecule`, `basis`,
  // `method`; `jk` must belong to `basis`. The references must outlive the
  // propagator.
  RealTimePropagator(const Molecule& molecule, const BasisSet& basis, const Method& method, const ScfResult& scf,
                     JKBuilder& jk);

  // Applies the pulse E(t) = kick * delta(t) (vector: strength times
  // direction, atomic units) to the ground state and propagates.
  RealTimeTrajectory propagate(const Vec3& kick, const RealTimeOptions& options = {}) const;

  // Fock matrix and total energy for a complex Hermitian density matrix
  // P = 2 C C^+ in the atomic orbital basis:
  //   F = h + J[Re P] + V_xc[Re P] - 1/2 (c_x K[P] + c_lr K_LR[P]),
  //   K[P] = K[Re P] + i K[Im P].
  // Public for tests.
  struct FockAndEnergy {
    Eigen::MatrixXcd fock;
    double energy = 0.0;
  };
  FockAndEnergy fock(const Eigen::MatrixXcd& density) const;

  Vec3 ground_state_dipole() const { return ground_state_dipole_; }
  Eigen::Index n_orbitals() const { return C_.cols(); }
  Eigen::Index n_occupied() const { return n_occ_; }

 private:
  const BasisSet& basis_;
  Functional functional_;
  JKBuilder& jk_;
  Eigen::MatrixXd C_;                       // ground-state orbitals: orthonormal basis of the propagation
  Eigen::Index n_occ_;
  Eigen::MatrixXd h_;                       // T + V
  std::array<Eigen::MatrixXd, 3> dipole_;   // atomic orbital basis
  Vec3 nuclear_dipole_{}, ground_state_dipole_{};
  double nuclear_repulsion_ = 0.0;
  std::unique_ptr<XcEvaluator> xc_;
  MolecularGrid grid_;
};

}  // namespace ellidft
