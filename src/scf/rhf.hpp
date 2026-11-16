#pragma once

#include <functional>
#include <vector>

#include <Eigen/Dense>

#include "basis/basis_set.hpp"
#include "core/molecule.hpp"
#include "grid/molecular_grid.hpp"
#include "scf/jk.hpp"
#include "xc/functional.hpp"

namespace ellidft {

enum class InitialGuess {
  AtomicDensities,   // superposition of spherically averaged atomic densities
  CoreHamiltonian,   // orbitals of the one-electron Hamiltonian
};

struct ScfIteration {
  int iteration = 0;
  double energy = 0.0;         // total energy, Hartree
  double energy_change = 0.0;  // to the previous iteration
  double error = 0.0;          // largest element of the commutator
};

struct ScfOptions {
  // Converged when the energy changes by less than energy_tolerance (Hartree)
  // between two iterations and the largest element of the commutator
  // F D S - S D F, in the orthonormal basis, is below error_tolerance.
  double energy_tolerance = 1e-10;
  double error_tolerance = 1e-8;
  int max_iterations = 200;
  std::size_t diis_vectors = 8;

  // Near-linear dependencies: eigenvectors of the overlap matrix with
  // eigenvalues below this threshold are removed from the orbital space
  // (canonical orthogonalisation). The reference code removes nothing by
  // default; for comparisons the threshold must be matched.
  double overlap_threshold = 1e-8;

  InitialGuess initial_guess = InitialGuess::AtomicDensities;

  // Start density from outside (restart from a checkpoint): if not empty it
  // replaces the initial guess. Must have the dimension of the basis.
  Eigen::MatrixXd start_density;

  // Called after every iteration with the density matrix that entered it
  // (for checkpoints and protocols). May be empty.
  std::function<void(const ScfIteration&, const Eigen::MatrixXd& density)> on_iteration;

  // With a J/K builder that screens with the density, J and K are updated
  // from the change of the density matrix; every full_rebuild_interval-th
  // iteration they are rebuilt from the full density.
  bool difference_densities = true;
  int full_rebuild_interval = 8;

  // Homogeneous external electric field F in atomic units. The electrons feel
  // the potential +F.r, the nuclei -Z F.R, so that E(F) = E(0) - mu.F - ...
  // with the dipole moment mu (nuclei minus electrons). Used for
  // finite-field checks of response properties.
  Vec3 electric_field{0.0, 0.0, 0.0};
};

struct ScfTimings {  // wall time in seconds
  double one_electron = 0.0;
  double initial_guess = 0.0;
  double coulomb_exchange = 0.0;
  double exchange_correlation = 0.0;
  double diagonalisation = 0.0;
  double total = 0.0;
};


struct ScfResult {
  bool converged = false;
  int iterations = 0;            // number of Fock matrix constructions
  int n_occupied = 0;
  int n_orbitals = 0;            // basis functions minus removed linear dependencies
  double energy_total = 0.0;     // electronic energy + nuclear repulsion, Hartree
  double energy_nuclear_repulsion = 0.0;  // plus the energy of the nuclei in an external field
  double energy_one_electron = 0.0;   // Tr[D (T + V)]
  double energy_two_electron = 0.0;   // 1/2 Tr[D (J - (c_x K + c_lr K_LR)/2)], exact-exchange part included
  double energy_xc = 0.0;             // density functional part, zero for Hartree-Fock
  double electrons_on_grid = 0.0;     // integral of the density on the grid (Kohn-Sham only)
  double energy_change = 0.0;    // last change of the total energy
  double error = 0.0;            // last commutator norm (largest element)
  Vec3 dipole{};                 // e * Bohr, nuclei minus electrons, origin (0,0,0)
  Eigen::VectorXd orbital_energies;
  Eigen::MatrixXd coefficients;  // columns are molecular orbitals
  Eigen::MatrixXd density;       // D = 2 C_occ C_occ^T
  Eigen::MatrixXd fock;          // F[D]
  Eigen::MatrixXd overlap;
  double overlap_min_eigenvalue = 0.0;
  ScfTimings timings;
  std::vector<ScfIteration> history;   // one entry per iteration, for protocols
};

// Superposition of atomic densities: for every atom a spherically averaged
// Hartree-Fock density in the atom's own basis functions (fractional
// occupation of the open shell), placed on the diagonal blocks of the
// molecular density matrix. Tr(D S) is the sum of the nuclear charges.
Eigen::MatrixXd atomic_density_guess(const Molecule& molecule, const BasisSet& basis);

// Method of a self-consistent field calculation: functional and grid.
struct Method {
  Functional functional = Functional::by_name("HF");
  GridSpec grid;  // used only if the functional has a density part
};

// Closed-shell Kohn-Sham or Hartree-Fock calculation:
//
//   F = T + V + J - (c_x K + c_lr K_LR) / 2 + V_xc,
//   E = Tr[D (T + V)] + 1/2 Tr[D (J - (c_x K + c_lr K_LR) / 2)] + E_xc + E_nuc,
//
// with the exact-exchange coefficients c_x, c_lr and the range separation
// parameter of the functional (K_LR: exchange with erf(omega r)/r).
//
// See ScfOptions for the start density, the treatment of linear
// dependencies and the use of density differences. Throws
// std::invalid_argument for an odd number of electrons or more occupied
// orbitals than independent basis functions.
ScfResult run_scf(const Molecule& molecule, const BasisSet& basis, const Method& method, JKBuilder& jk,
                  const ScfOptions& options = {});
ScfResult run_scf(const Molecule& molecule, const BasisSet& basis, const Method& method,
                  const ScfOptions& options = {});

// Restricted Hartree-Fock for a closed-shell molecule (Hernández Vera 2025,
// Eqs. 21-24 and section 3.2): run_scf with exact exchange only.
ScfResult run_rhf(const Molecule& molecule, const BasisSet& basis, JKBuilder& jk,
                  const ScfOptions& options = {});

// Same with the integral-direct four-centre J/K builder.
ScfResult run_rhf(const Molecule& molecule, const BasisSet& basis, const ScfOptions& options = {});

}  // namespace ellidft
