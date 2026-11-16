#include "realtime/propagator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "integrals/one_electron.hpp"
#include "scf/xc_matrix.hpp"

namespace ellidft {

namespace {

using Complex = std::complex<double>;

// exp(-i * factor * H) for a Hermitian matrix H through its eigenvectors.
// The result is unitary to rounding for any step length.
Eigen::MatrixXcd exponential(const Eigen::MatrixXcd& hermitian, double factor) {
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(hermitian);
  if (eig.info() != Eigen::Success) throw std::runtime_error("real-time propagation: diagonalisation failed");
  const Eigen::VectorXcd phases =
      (Complex(0.0, -factor) * eig.eigenvalues().cast<Complex>()).array().exp().matrix();
  return eig.eigenvectors() * phases.asDiagonal() * eig.eigenvectors().adjoint();
}

}  // namespace

RealTimePropagator::RealTimePropagator(const Molecule& molecule, const BasisSet& basis, const Method& method,
                                       const ScfResult& scf, JKBuilder& jk)
    : basis_(basis), functional_(method.functional), jk_(jk), C_(scf.coefficients), n_occ_(scf.n_occupied) {
  if (!scf.converged) throw std::invalid_argument("RealTimePropagator: the SCF is not converged");
  if (n_occ_ <= 0 || n_occ_ > C_.cols()) throw std::invalid_argument("RealTimePropagator: no occupied orbitals");
  h_ = kinetic_matrix(basis) + nuclear_attraction_matrix(basis, molecule);
  dipole_ = dipole_matrices(basis);
  nuclear_repulsion_ = molecule.nuclear_repulsion_energy();
  for (const Atom& atom : molecule.atoms)
    for (std::size_t k = 0; k < 3; ++k) nuclear_dipole_[k] += atom.Z * atom.r[k];
  for (std::size_t k = 0; k < 3; ++k)
    ground_state_dipole_[k] = nuclear_dipole_[k] - scf.density.cwiseProduct(dipole_[k]).sum();
  if (functional_.has_density_part()) {
    xc_ = make_libxc_evaluator(functional_);
    grid_ = build_molecular_grid(molecule, method.grid);
  }
}

RealTimePropagator::FockAndEnergy RealTimePropagator::fock(const Eigen::MatrixXcd& density) const {
  // With complex orbitals the exchange matrix is
  //   K_mu,nu = sum_{lambda sigma} (mu lambda | nu sigma) P_lambda,sigma,  P = 2 C C^+,
  // (from K chi_nu = sum_i psi_i <psi_i| 1/r12 |chi_nu> with real basis
  // functions), i.e. the K of the J/K interface applied to the real and the
  // imaginary part of P. Re P is symmetric, Im P antisymmetric; J[Im P]
  // vanishes and K[Im P] is antisymmetric, so F is Hermitian. The density
  // rho = sum P_mu,nu chi_mu chi_nu contains only Re P.
  const double c_x = functional_.exact_exchange, c_lr = functional_.long_range_exchange;
  const Eigen::MatrixXd re = density.real(), im = density.imag();
  const bool exchange = c_x != 0.0 || c_lr != 0.0;

  std::vector<Eigen::MatrixXd> J, K, K_lr;
  if (c_x != 0.0) jk_.build({re, im}, &J, &K);
  else jk_.build({re}, &J, nullptr);
  if (c_lr != 0.0) jk_.build_long_range_exchange({re, im}, functional_.omega, K_lr);

  FockAndEnergy result;
  Eigen::MatrixXd f_re = h_ + J[0];
  Eigen::MatrixXd f_im = Eigen::MatrixXd::Zero(re.rows(), re.cols());
  if (c_x != 0.0) {
    f_re -= 0.5 * c_x * K[0];
    f_im -= 0.5 * c_x * K[1];
  }
  if (c_lr != 0.0) {
    f_re -= 0.5 * c_lr * K_lr[0];
    f_im -= 0.5 * c_lr * K_lr[1];
  }
  // E = Tr[P h] + 1/2 Tr[P (J - 1/2 (c_x K + c_lr K_LR))] + E_xc + E_nuc;
  // Tr[P K] = sum P_nu,mu K_mu,nu = sum Re P * K[Re P] + sum Im P * K[Im P]
  // (both sums element-wise; the mixed terms vanish by symmetry).
  double energy = re.cwiseProduct(h_).sum() + 0.5 * re.cwiseProduct(J[0]).sum() + nuclear_repulsion_;
  if (exchange) {
    const Eigen::MatrixXd g_re = f_re - h_ - J[0];   // -1/2 (c_x K + c_lr K_LR)[Re P]
    energy += 0.5 * (re.cwiseProduct(g_re).sum() + im.cwiseProduct(f_im).sum());
  }
  if (xc_) {
    const XcIntegrals xc = integrate_xc(basis_, grid_, *xc_, re);
    f_re += xc.matrix;
    energy += xc.energy;
  }
  result.fock.resize(re.rows(), re.cols());
  result.fock.real() = f_re;
  result.fock.imag() = f_im;
  result.energy = energy;
  return result;
}

RealTimeTrajectory RealTimePropagator::propagate(const Vec3& kick, const RealTimeOptions& options) const {
  if (!(options.time_step > 0.0)) throw std::invalid_argument("real-time propagation: time step must be positive");
  if (options.n_steps < 1) throw std::invalid_argument("real-time propagation: needs at least one step");
  if (options.max_fock_builds_per_step < 1)
    throw std::invalid_argument("real-time propagation: needs at least one Fock build per step");

  const Eigen::Index n = C_.cols();
  const Eigen::MatrixXcd C_ao = C_.cast<Complex>();
  const Complex i_unit(0.0, 1.0);
  const double dt = options.time_step;
  // Occupied orbitals as columns in the basis of the ground-state orbitals.
  Eigen::MatrixXcd orbitals = Eigen::MatrixXcd::Zero(n, n_occ_);
  for (Eigen::Index i = 0; i < n_occ_; ++i) orbitals(i, i) = 1.0;

  std::array<Eigen::MatrixXd, 3> dipole_mo;
  for (std::size_t k = 0; k < 3; ++k) dipole_mo[k] = C_.transpose() * dipole_[k] * C_;

  // Pulse: H(t) = H0 - mu.u F(t) with F(t) = kappa delta(t) (Hauge et al.
  // 2023, Eqs. 3, 4, 21). Integrating the Schroedinger equation across the
  // pulse gives psi(0+) = exp(+i kappa mu.u) psi(0-); for electrons mu = -r,
  // so every orbital is multiplied by exp(-i kappa u.r), here with the
  // dipole operator projected onto the orbital space.
  Eigen::MatrixXd kick_operator = Eigen::MatrixXd::Zero(n, n);
  for (std::size_t k = 0; k < 3; ++k)
    if (kick[k] != 0.0) kick_operator += kick[k] * dipole_mo[k];
  if (kick_operator.cwiseAbs().maxCoeff() > 0.0) orbitals = exponential(kick_operator.cast<Complex>(), 1.0) * orbitals;

  // Field during the propagation: the electrons feel +F(t).r.
  auto field_term = [&](double time) {
    Eigen::MatrixXcd V = Eigen::MatrixXcd::Zero(n, n);
    if (options.field) {
      const Vec3 f = options.field(time);
      for (std::size_t k = 0; k < 3; ++k)
        if (f[k] != 0.0) V += f[k] * dipole_mo[k].cast<Complex>();
    }
    return V;
  };

  struct State {
    Eigen::MatrixXcd fock;   // orthonormal basis, without the field
    double energy = 0.0;
    Vec3 dipole{};
  };
  RealTimeTrajectory trajectory;
  trajectory.time_step = dt;
  auto evaluate = [&](const Eigen::MatrixXcd& occupied) {
    const Eigen::MatrixXcd c = C_ao * occupied;
    const Eigen::MatrixXcd density = 2.0 * c * c.adjoint();
    const FockAndEnergy fe = fock(density);
    ++trajectory.fock_builds;
    State state;
    const Eigen::MatrixXcd f = C_ao.adjoint() * fe.fock * C_ao;
    state.fock = 0.5 * (f + f.adjoint());
    state.energy = fe.energy;
    for (std::size_t k = 0; k < 3; ++k)
      state.dipole[k] = nuclear_dipole_[k] - density.real().cwiseProduct(dipole_[k]).sum();
    return state;
  };

  State current = evaluate(orbitals);
  trajectory.dipole.push_back(current.dipole);
  trajectory.energy.push_back(current.energy);
  std::vector<Eigen::MatrixXcd> history;   // Fock matrices at t_n-1, t_n-2 (in this order)

  // Weights of the polynomial through values at the times `nodes` (in units
  // of the step, relative to t_n), evaluated at `at` (Lagrange form).
  auto lagrange = [](const std::vector<double>& nodes, double at) {
    std::vector<double> weights(nodes.size(), 1.0);
    for (std::size_t a = 0; a < nodes.size(); ++a)
      for (std::size_t b = 0; b < nodes.size(); ++b)
        if (a != b) weights[a] *= (at - nodes[b]) / (nodes[a] - nodes[b]);
    return weights;
  };
  auto note = [&](bool converged, int builds) {
    trajectory.steps_converged = trajectory.steps_converged && converged;
    trajectory.max_fock_builds_in_a_step = std::max(trajectory.max_fock_builds_in_a_step, builds);
  };

  // --- Midpoint step ------------------------------------------------------
  // From time t_a over the length s (may be negative):
  //   C_b = exp(-i s F_M) C_a,   F_M = (F[C_a] + F[C_b]) / 2 + V(t_a + s/2)
  // (Herbert 2022, Eqs. 78-80, there for the density matrix), made
  // self-consistent in F[C_b]. The step is symmetric in its two ends. For
  // Hartree-Fock it conserves the energy exactly when V = 0: the energy is
  // quadratic in the density matrix, hence E_b - E_a = Tr[F_M (P_b - P_a)],
  // which vanishes because U = exp(-i s F_M) commutes with F_M.
  auto midpoint_step = [&](double t_a, double s, const Eigen::MatrixXcd& from, const State& a, Eigen::MatrixXcd trial,
                           Eigen::MatrixXcd& to, State& b) {
    const Eigen::MatrixXcd V = field_term(t_a + 0.5 * s);
    bool converged = false;
    int builds = 0;
    while (builds < options.max_fock_builds_per_step) {
      const Eigen::MatrixXcd H = s * (0.5 * (a.fock + trial) + V);
      to = exponential(0.5 * (H + H.adjoint()), 1.0) * from;
      b = evaluate(to);
      ++builds;
      const double change = (b.fock - trial).cwiseAbs().maxCoeff();
      trial = b.fock;
      if (change <= options.fock_tolerance) {
        converged = true;
        break;
      }
    }
    note(converged, builds);
  };

  // --- Magnus step --------------------------------------------------------
  // For Y' = A(t) Y the fourth-order Magnus integrator is
  //   Y_n+1 = exp(Omega) Y_n,  Omega = (h/2) (A1 + A2) - (sqrt(3) h^2 / 12) [A1, A2],
  //   A1,2 = A(t_n + (1/2 -+ sqrt(3)/6) h)
  // (Blanes et al. 2008, Eqs. 252-254). With A = -i F this is
  //   U = exp(-i H),  H = (h/2) (F1 + F2) + i (sqrt(3) h^2 / 12) [F1, F2],
  // and H is Hermitian.
  //
  // F inside the step depends on the orbitals at its end. This is handled
  // here, not in the sources: the field-free part of F(t) is represented by
  // the polynomial through the Fock matrices at t_n+1 (trial), t_n and up to
  // two earlier times (cubic, error of fourth order), and the step is
  // repeated with the Fock matrix of the orbitals it produces until that no
  // longer changes. The first trial is the extrapolation of the earlier
  // matrices. The external field is added at the Gauss points exactly.
  //
  // Directly after the pulse there are no earlier times, and a polynomial of
  // lower degree would leave an error of third order in the whole
  // trajectory. Therefore the first three steps are refined together: with
  // the cubic through the Fock matrices at t_0 ... t_3 the three steps are
  // repeated until these matrices no longer change.
  const double c1 = 0.5 - std::sqrt(3.0) / 6.0, c2 = 0.5 + std::sqrt(3.0) / 6.0;
  // H for the step that starts at `t_start`; F given at `nodes` (in units of
  // the step, relative to t_start).
  auto magnus_generator = [&](double t_start, const std::vector<double>& nodes,
                              const std::vector<const Eigen::MatrixXcd*>& matrices) {
    auto at = [&](double c) {
      const std::vector<double> weights = lagrange(nodes, c);
      Eigen::MatrixXcd sum = field_term(t_start + c * dt);
      for (std::size_t k = 0; k < nodes.size(); ++k) sum += weights[k] * *matrices[k];
      return sum;
    };
    const Eigen::MatrixXcd F1 = at(c1), F2 = at(c2);
    const Eigen::MatrixXcd H = 0.5 * dt * (F1 + F2) + i_unit * (std::sqrt(3.0) * dt * dt / 12.0) * (F1 * F2 - F2 * F1);
    return Eigen::MatrixXcd(0.5 * (H + H.adjoint()));
  };

  // --- Composition --------------------------------------------------------
  // "Triple jump": from a symmetric scheme of order 2 a scheme of order 4
  // follows by three steps of lengths w1 h, w0 h, w1 h (Wikipedia,
  // "Symplectic integrator", method of triple jump; coefficients
  // w1 = 1 / (2 - 2^(1/3)), w0 = -2^(1/3) / (2 - 2^(1/3)) from Wikipedia,
  // "Leapfrog integration", Yoshida algorithm). Applied here to the midpoint
  // step above; the middle step runs backwards in time.
  const double w1 = 1.0 / (2.0 - std::cbrt(2.0)), w0 = -std::cbrt(2.0) / (2.0 - std::cbrt(2.0));

  const bool refine_start = options.propagator == Propagator::Magnus4 && options.n_steps >= 3;
  std::vector<Eigen::MatrixXcd> start_orbitals = {orbitals};
  std::vector<State> start_states = {current};
  // previous point of the trajectory, for the first trial of a midpoint step
  bool have_previous = false;
  double previous_time = 0.0;
  Eigen::MatrixXcd previous_fock;

  for (int step = 0; step < options.n_steps; ++step) {
    const double t_n = dt * step;
    State next;
    Eigen::MatrixXcd next_orbitals;

    if (options.propagator == Propagator::Composition4) {
      double t_a = t_n;
      Eigen::MatrixXcd from = orbitals;
      State a = current;
      for (double s : {w1 * dt, w0 * dt, w1 * dt}) {
        // first trial: straight line through the last two points in time
        Eigen::MatrixXcd trial = a.fock;
        if (have_previous) trial += (s / (t_a - previous_time)) * (a.fock - previous_fock);
        midpoint_step(t_a, s, from, a, trial, next_orbitals, next);
        have_previous = true;
        previous_time = t_a;
        previous_fock = a.fock;
        t_a += s;
        from = next_orbitals;
        a = next;
      }
    } else {
      // known Fock matrices: at node time 0 (current), -1, -2 (history)
      std::vector<double> nodes = {0.0};
      std::vector<const Eigen::MatrixXcd*> matrices = {&current.fock};
      for (std::size_t k = 0; k < history.size(); ++k) {
        nodes.push_back(-static_cast<double>(k + 1));
        matrices.push_back(&history[k]);
      }
      // first trial for F at the end of the step: extrapolation
      Eigen::MatrixXcd trial = Eigen::MatrixXcd::Zero(n, n);
      {
        const std::vector<double> weights = lagrange(nodes, 1.0);
        for (std::size_t k = 0; k < nodes.size(); ++k) trial += weights[k] * *matrices[k];
      }
      if (options.propagator == Propagator::Midpoint) {
        midpoint_step(t_n, dt, orbitals, current, trial, next_orbitals, next);
      } else {
        nodes.push_back(1.0);
        matrices.push_back(&trial);
        bool converged = false;
        int builds = 0;
        while (builds < options.max_fock_builds_per_step) {
          next_orbitals = exponential(magnus_generator(t_n, nodes, matrices), 1.0) * orbitals;
          next = evaluate(next_orbitals);
          ++builds;
          const double change = (next.fock - trial).cwiseAbs().maxCoeff();
          trial = next.fock;
          if (change <= options.fock_tolerance) {
            converged = true;
            break;
          }
        }
        note(converged, builds);
      }
      history.insert(history.begin(), current.fock);
      if (history.size() > 2) history.pop_back();
    }

    orbitals = next_orbitals;
    current = next;
    trajectory.dipole.push_back(current.dipole);
    trajectory.energy.push_back(current.energy);

    if (refine_start && step < 3) {
      start_orbitals.push_back(orbitals);
      start_states.push_back(current);
    }
    if (refine_start && step == 2) {
      bool refined = false;
      for (int pass = 0; pass < options.max_fock_builds_per_step && !refined; ++pass) {
        double change = 0.0;
        for (int k = 0; k < 3; ++k) {
          const std::vector<double> block_nodes = {0.0 - k, 1.0 - k, 2.0 - k, 3.0 - k};
          const std::vector<const Eigen::MatrixXcd*> block = {&start_states[0].fock, &start_states[1].fock,
                                                              &start_states[2].fock, &start_states[3].fock};
          const auto to = static_cast<std::size_t>(k) + 1;
          start_orbitals[to] = exponential(magnus_generator(dt * k, block_nodes, block), 1.0) * start_orbitals[to - 1];
          State improved = evaluate(start_orbitals[to]);
          change = std::max(change, (improved.fock - start_states[to].fock).cwiseAbs().maxCoeff());
          start_states[to] = std::move(improved);
        }
        refined = change <= options.fock_tolerance;
      }
      trajectory.steps_converged = trajectory.steps_converged && refined;
      for (std::size_t k = 1; k <= 3; ++k) {
        trajectory.dipole[k] = start_states[k].dipole;
        trajectory.energy[k] = start_states[k].energy;
      }
      orbitals = start_orbitals[3];
      current = start_states[3];
      history = {start_states[2].fock, start_states[1].fock};
    }
  }
  for (double energy : trajectory.energy)
    trajectory.energy_drift = std::max(trajectory.energy_drift, std::abs(energy - trajectory.energy.front()));

  trajectory.orthonormality_error =
      (orbitals.adjoint() * orbitals - Eigen::MatrixXcd::Identity(n_occ_, n_occ_)).cwiseAbs().maxCoeff();
  return trajectory;
}

}  // namespace ellidft
