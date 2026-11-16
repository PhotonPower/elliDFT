#include "scf/rhf.hpp"

#include <chrono>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "integrals/one_electron.hpp"
#include "scf/diis.hpp"
#include "scf/xc_matrix.hpp"
#include "xc/xc_evaluator.hpp"

namespace ellidft {

namespace {

class Stopwatch {
 public:
  Stopwatch() : start_(std::chrono::steady_clock::now()) {}
  double seconds() const {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start_).count();
  }

 private:
  std::chrono::steady_clock::time_point start_;
};

// How the orbitals are occupied.
enum class Occupation {
  ClosedShell,          // the lowest N/2 orbitals doubly
  FractionalDegenerate  // atoms: electrons spread evenly over a degenerate level
};

// Occupation numbers (0..2) for orbital energies in ascending order.
Eigen::VectorXd occupation_numbers(const Eigen::VectorXd& energies, int n_electrons, Occupation occupation) {
  Eigen::VectorXd occ = Eigen::VectorXd::Zero(energies.size());
  if (occupation == Occupation::ClosedShell) {
    occ.head(n_electrons / 2).setConstant(2.0);
    return occ;
  }
  // Levels closer than this are treated as one degenerate level. In a
  // spherical atom the members of a shell are degenerate to rounding.
  const double degeneracy_tolerance = 1e-6;
  double remaining = n_electrons;
  Eigen::Index i = 0;
  while (remaining > 0.0 && i < energies.size()) {
    Eigen::Index j = i;
    while (j + 1 < energies.size() && energies(j + 1) - energies(i) < degeneracy_tolerance) ++j;
    const auto size = static_cast<double>(j - i + 1);
    const double electrons = std::min(remaining, 2.0 * size);
    for (Eigen::Index k = i; k <= j; ++k) occ(k) = electrons / size;
    remaining -= electrons;
    i = j + 1;
  }
  if (remaining > 0.0) throw std::invalid_argument("SCF: more electrons than the orbitals can hold");
  return occ;
}

ScfResult solve(const Molecule& molecule, const BasisSet& basis, const Method& method, JKBuilder& jk,
                const ScfOptions& options, Occupation occupation) {
  const Stopwatch total_time;
  const int n_electrons = molecule.n_electrons();
  if (n_electrons < 0) throw std::invalid_argument("run_scf: negative number of electrons");
  if (occupation == Occupation::ClosedShell && n_electrons % 2 != 0)
    throw std::invalid_argument("run_scf: closed-shell calculation needs an even number of electrons, got " +
                                std::to_string(n_electrons));
  const Eigen::Index n = basis.n_functions();

  const Functional& functional = method.functional;
  const double c_x = functional.exact_exchange;
  const double c_lr = functional.long_range_exchange;
  const bool use_exchange = c_x != 0.0;
  const bool use_long_range = c_lr != 0.0;
  if (use_long_range && !(functional.omega > 0.0))
    throw std::invalid_argument("run_scf: long-range exchange needs a positive range separation parameter");
  std::unique_ptr<XcEvaluator> xc;
  MolecularGrid grid;
  if (functional.has_density_part()) {
    xc = make_libxc_evaluator(functional);
    grid = build_molecular_grid(molecule, method.grid);
  }

  ScfResult result;
  result.n_occupied = n_electrons / 2;
  result.energy_nuclear_repulsion = molecule.nuclear_repulsion_energy();

  const Stopwatch one_electron_time;
  const Eigen::MatrixXd S = overlap_matrix(basis);
  Eigen::MatrixXd H = kinetic_matrix(basis) + nuclear_attraction_matrix(basis, molecule);
  const Vec3& field = options.electric_field;
  if (field[0] != 0.0 || field[1] != 0.0 || field[2] != 0.0) {
    // Interaction with a homogeneous field F: -mu.F with mu = sum Z R - r,
    // i.e. +F.r for the electrons and -Z F.R for the nuclei.
    const auto dipole = dipole_matrices(basis);
    for (std::size_t k = 0; k < 3; ++k) {
      H += field[k] * dipole[k];
      for (const Atom& atom : molecule.atoms) result.energy_nuclear_repulsion -= atom.Z * atom.r[k] * field[k];
    }
  }
  result.timings.one_electron = one_electron_time.seconds();

  // Canonical orthogonalisation X = U s^(-1/2), restricted to the eigenvectors
  // of S with eigenvalues above the threshold; X^T S X = 1. Without removed
  // vectors this spans the same space as the symmetric orthogonalisation
  // U s^(-1/2) U^T of Hernández Vera 2025, section 3.2, step 3.
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> overlap_eigen(S);
  result.overlap_min_eigenvalue = overlap_eigen.eigenvalues().minCoeff();
  Eigen::Index n_removed = 0;
  while (n_removed < n && overlap_eigen.eigenvalues()(n_removed) < options.overlap_threshold) ++n_removed;
  const Eigen::Index m = n - n_removed;  // eigenvalues are in ascending order
  result.n_orbitals = static_cast<int>(m);
  if (occupation == Occupation::ClosedShell && n_electrons / 2 > m)
    throw std::invalid_argument("run_scf: more occupied orbitals than independent basis functions");
  const Eigen::MatrixXd X = overlap_eigen.eigenvectors().rightCols(m) *
                            overlap_eigen.eigenvalues().tail(m).cwiseInverse().cwiseSqrt().asDiagonal();

  // Orbitals and density from a Fock matrix: F' = X^T F X, F' C' = C' eps,
  // C = X C', D = sum_i n_i C_i C_i^T with n_i = 2 for occupied orbitals
  // (section 3.2, steps 4c, 4d; Eq. 24).
  Eigen::VectorXd orbital_energies;
  Eigen::MatrixXd C;
  auto diagonalize = [&](const Eigen::MatrixXd& F) {
    const Stopwatch time;
    const Eigen::MatrixXd F_orthogonal = X.transpose() * F * X;
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(F_orthogonal);
    orbital_energies = eigen.eigenvalues();
    C = X * eigen.eigenvectors();
    result.timings.diagonalisation += time.seconds();
  };
  auto density_from_orbitals = [&]() -> Eigen::MatrixXd {
    const Eigen::VectorXd occ = occupation_numbers(orbital_energies, n_electrons, occupation);
    return C * occ.asDiagonal() * C.transpose();
  };

  // Start density.
  const Stopwatch guess_time;
  Eigen::MatrixXd D;
  if (options.start_density.size() > 0) {
    if (options.start_density.rows() != n || options.start_density.cols() != n)
      throw std::invalid_argument("run_scf: the start density does not have the dimension of the basis");
    D = options.start_density;
  } else if (options.initial_guess == InitialGuess::AtomicDensities && occupation == Occupation::ClosedShell &&
      molecule.n_atoms() > 0) {
    D = atomic_density_guess(molecule, basis);
  } else {
    diagonalize(H);
    D = density_from_orbitals();
  }
  result.timings.initial_guess = guess_time.seconds();

  Diis diis(options.diis_vectors);
  const bool incremental = options.difference_densities && jk.screens_with_density();
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(n, n), K = Eigen::MatrixXd::Zero(n, n), D_previous;
  Eigen::MatrixXd K_lr = Eigen::MatrixXd::Zero(n, n);
  std::vector<Eigen::MatrixXd> J_part, K_part, K_lr_part;
  double previous_energy = 0.0;

  for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
    // Coulomb and exchange matrices are linear in the density. With a
    // screening builder only the change of the density is processed, which
    // lets more and more integrals drop out as the SCF converges.
    const Stopwatch jk_time;
    const bool full_build = !incremental || iteration == 1 ||
                            (options.full_rebuild_interval > 0 && (iteration - 1) % options.full_rebuild_interval == 0);
    if (full_build) {
      jk.build({D}, &J_part, use_exchange ? &K_part : nullptr);
      J = J_part[0];
      if (use_exchange) K = K_part[0];
      if (use_long_range) {
        jk.build_long_range_exchange({D}, functional.omega, K_lr_part);
        K_lr = K_lr_part[0];
      }
    } else {
      const Eigen::MatrixXd dD = D - D_previous;
      jk.build({dD}, &J_part, use_exchange ? &K_part : nullptr);
      J += J_part[0];
      if (use_exchange) K += K_part[0];
      if (use_long_range) {
        jk.build_long_range_exchange({dD}, functional.omega, K_lr_part);
        K_lr += K_lr_part[0];
      }
    }
    D_previous = D;
    result.timings.coulomb_exchange += jk_time.seconds();

    // Hartree-Fock: F = H + J - K/2 (Eq. 23), E = 1/2 Tr[D (H + F)]
    // (section 3.2, step 4e). In Kohn-Sham theory the exchange term is
    // scaled by the fraction c_x of exact exchange and the density functional
    // part adds V_xc to F and E_xc to E.
    // Range-separated hybrids add the long-range exchange K_LR (operator
    // erf(omega r)/r) with its own coefficient: G = J - (c_x K + c_lr K_LR)/2.
    Eigen::MatrixXd G = J;
    if (use_exchange) G -= 0.5 * c_x * K;
    if (use_long_range) G -= 0.5 * c_lr * K_lr;
    Eigen::MatrixXd F = H + G;
    const double energy_one = D.cwiseProduct(H).sum();
    const double energy_two = 0.5 * D.cwiseProduct(G).sum();
    double energy_xc = 0.0;
    if (xc) {
      const Stopwatch xc_time;
      const XcIntegrals xc_integrals = integrate_xc(basis, grid, *xc, D);
      F += xc_integrals.matrix;
      energy_xc = xc_integrals.energy;
      result.electrons_on_grid = xc_integrals.n_electrons;
      result.timings.exchange_correlation += xc_time.seconds();
    }
    const double energy = energy_one + energy_two + energy_xc + result.energy_nuclear_repulsion;

    // Commutator error (Eq. 65), in the orthonormal basis.
    const Eigen::MatrixXd FDS = F * D * S;
    const Eigen::MatrixXd error = X.transpose() * (FDS - FDS.transpose()) * X;

    result.iterations = iteration;
    result.energy_total = energy;
    result.energy_one_electron = energy_one;
    result.energy_two_electron = energy_two;
    result.energy_xc = energy_xc;
    result.energy_change = energy - previous_energy;
    result.error = error.cwiseAbs().maxCoeff();
    result.history.push_back({iteration, energy, result.energy_change, result.error});
    if (options.on_iteration) options.on_iteration(result.history.back(), D);
    result.density = D;
    result.fock = F;

    if (iteration > 1 && std::abs(result.energy_change) < options.energy_tolerance &&
        result.error < options.error_tolerance) {
      result.converged = true;
      break;
    }
    previous_energy = energy;
    if (iteration == options.max_iterations) break;

    diagonalize(diis.extrapolate(F, error));
    D = density_from_orbitals();
  }

  // Orbitals of the final Fock matrix F[D]. The density and the energy are
  // those of the last iteration; at convergence the two are consistent.
  diagonalize(result.fock);
  result.orbital_energies = orbital_energies;
  result.coefficients = C;
  result.overlap = S;

  // Dipole moment: point nuclei minus the electron density.
  const auto dipole = dipole_matrices(basis);
  for (std::size_t k = 0; k < 3; ++k) {
    double value = 0.0;
    for (const Atom& atom : molecule.atoms) value += atom.Z * atom.r[k];
    result.dipole[k] = value - result.density.cwiseProduct(dipole[k]).sum();
  }
  result.timings.total = total_time.seconds();
  return result;
}

bool same_shell_data(const Shell& a, const Shell& b) {
  return a.l == b.l && a.spherical == b.spherical && a.exponents == b.exponents && a.coefficients == b.coefficients;
}

}  // namespace

Eigen::MatrixXd atomic_density_guess(const Molecule& molecule, const BasisSet& basis) {
  const Eigen::Index n = basis.n_functions();
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n, n);
  const auto& shells = basis.shells();

  // Shells and first basis function of every atom (the shells of an atom are
  // contiguous in the basis).
  struct AtomBlock {
    std::vector<Shell> shells;
    Eigen::Index first = 0, size = 0;
  };
  std::vector<AtomBlock> blocks(molecule.n_atoms());
  for (std::size_t s = 0; s < shells.size(); ++s) {
    AtomBlock& block = blocks.at(static_cast<std::size_t>(shells[s].atom));
    if (block.shells.empty()) block.first = basis.first_function(s);
    block.shells.push_back(shells[s]);
    block.size += shells[s].n_functions();
  }

  // One atomic calculation per kind of atom (element and basis functions).
  struct Solved {
    int Z;
    std::vector<Shell> shells;
    Eigen::MatrixXd density;
  };
  std::vector<Solved> solved;
  ScfOptions atom_options;
  atom_options.initial_guess = InitialGuess::CoreHamiltonian;
  atom_options.energy_tolerance = 1e-8;
  atom_options.error_tolerance = 1e-6;

  for (std::size_t a = 0; a < molecule.n_atoms(); ++a) {
    const AtomBlock& block = blocks[a];
    if (block.shells.empty()) continue;  // atom without basis functions
    const Solved* match = nullptr;
    for (const Solved& candidate : solved) {
      if (candidate.Z != molecule.atoms[a].Z || candidate.shells.size() != block.shells.size()) continue;
      bool equal = true;
      for (std::size_t s = 0; s < block.shells.size() && equal; ++s)
        equal = same_shell_data(candidate.shells[s], block.shells[s]);
      if (equal) {
        match = &candidate;
        break;
      }
    }
    if (match == nullptr) {
      Molecule atom;
      atom.atoms = {molecule.atoms[a]};
      std::vector<Shell> atom_shells = block.shells;
      for (Shell& shell : atom_shells) shell.atom = 0;
      const BasisSet atom_basis = BasisSet::from_shells(basis.name(), basis.spherical(), atom_shells);
      DirectJK atom_jk(atom_basis);
      // Spherically averaged Hartree-Fock: the electrons of the open shell are
      // spread evenly over its degenerate orbitals, and the closed-shell
      // energy expression is used with these occupation numbers. This is a
      // model for a start density, not a physical state of the atom.
      const ScfResult atom_scf = solve(atom, atom_basis, Method{}, atom_jk, atom_options, Occupation::FractionalDegenerate);
      solved.push_back({molecule.atoms[a].Z, block.shells, atom_scf.density});
      match = &solved.back();
    }
    D.block(block.first, block.first, block.size, block.size) = match->density;
  }
  return D;
}

ScfResult run_rhf(const Molecule& molecule, const BasisSet& basis, const ScfOptions& options) {
  return run_scf(molecule, basis, Method{}, options);
}

ScfResult run_rhf(const Molecule& molecule, const BasisSet& basis, JKBuilder& jk, const ScfOptions& options) {
  return run_scf(molecule, basis, Method{}, jk, options);
}

ScfResult run_scf(const Molecule& molecule, const BasisSet& basis, const Method& method,
                  const ScfOptions& options) {
  DirectJK jk(basis);
  return run_scf(molecule, basis, method, jk, options);
}

ScfResult run_scf(const Molecule& molecule, const BasisSet& basis, const Method& method, JKBuilder& jk,
                  const ScfOptions& options) {
  return solve(molecule, basis, method, jk, options, Occupation::ClosedShell);
}

}  // namespace ellidft
