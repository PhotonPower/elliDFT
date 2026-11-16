// ellidft: the main program. Reads a control file (JSON), runs the ground
// state and the linear response, and writes the results as JSON in the layout
// of the specification, section 8, together with a protocol (settings,
// versions, iterations, times).
//
// Usage: ellidft <control.json>
//
//   {
//     "name": "styrene_4mer",                 // base name of the output files (default: name of the control file)
//     "geometry": "styrene_4mer.xyz", "bohr": false, "charge": 0,
//     "basis": "def2-SVPD", "aux_basis": "def2-universal-JKFIT",
//     "functional": "CAM-B3LYP", "grid": [75, 302],
//     "scf": {"energy_tolerance": 1e-10, "error_tolerance": 1e-8, "max_iterations": 200, "guess": "atomic"},
//     "response": {"wavelengths_nm": [589, 633], "frequencies_hartree": [], "gamma_hartree": 0.0,
//                  "states": 0, "tda": false},
//     "threads": 0,                           // 0: all
//     "output_directory": ".",
//     "restart": true                         // continue from <name>.chk if it fits (default)
//   }
//
// Only "geometry" and "basis" are required. File names are relative to the
// control file. Writes <name>.json and <name>.log; the protocol also goes to
// standard error. The JSON file is what ellidft-optics reads.
//
// Checkpoint: <name>.chk holds the density matrix (rewritten after every SCF
// iteration) and every finished polarisability. A later run with the same
// molecule, basis and method starts the SCF from that density and does not
// repeat finished frequencies. Excitation energies are not stored.

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/elements.hpp"
#include "core/parallel.hpp"
#include "core/units.hpp"
#include "core/xyz.hpp"
#include "df/df_jk.hpp"
#include "optics/optics.hpp"
#include "response/excitations.hpp"
#include "response/polarizability.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"
#include "xc/xc_evaluator.hpp"

using namespace ellidft;
using nlohmann::json;
using nlohmann::ordered_json;

namespace {

void check_keys(const json& node, const std::string& where, std::initializer_list<const char*> allowed) {
  if (!node.is_object()) throw std::runtime_error(where + ": expected an object");
  for (const auto& [key, unused] : node.items()) {
    bool known = key == "comment";
    for (const char* a : allowed) known = known || key == a;
    if (!known) throw std::runtime_error(where + ": unknown key \"" + key + "\"");
  }
}

ordered_json matrix_json(const Eigen::Matrix3d& m) {
  ordered_json rows = ordered_json::array();
  for (int k = 0; k < 3; ++k) rows.push_back({m(k, 0), m(k, 1), m(k, 2)});
  return rows;
}

// Writes every line to the protocol file and to standard error.
class Protocol {
 public:
  explicit Protocol(const std::filesystem::path& path) : file_(path) {
    if (!file_) throw std::runtime_error("cannot write " + path.string());
  }
  void line(const std::string& text = "") {
    file_ << text << "\n";
    file_.flush();
    std::cerr << text << "\n";
  }

 private:
  std::ofstream file_;
};

std::string now() {
  const std::time_t t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char buffer[32];
  std::strftime(buffer, sizeof buffer, "%Y-%m-%d %H:%M:%S", std::localtime(&t));
  return buffer;
}

double seconds_since(const std::chrono::steady_clock::time_point& start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

template <typename... Args>
std::string format(const char* pattern, Args... args) {
  char buffer[512];
  std::snprintf(buffer, sizeof buffer, pattern, args...);
  return buffer;
}

// Checkpoint file: one line of JSON (what the data belong to, state of the
// SCF, finished response results), then the density matrix as raw doubles.
struct Checkpoint {
  std::string fingerprint;
  bool scf_converged = false;
  json response = json::object();   // key -> {"re": 3x3, "im": 3x3, "iterations": n}
  Eigen::MatrixXd density;
};

void write_checkpoint(const std::filesystem::path& path, const Checkpoint& checkpoint) {
  // written under another name first, so that an interrupted write leaves the old file intact
  std::filesystem::path temporary = path;
  temporary += ".tmp";
  {
    std::ofstream file(temporary, std::ios::binary);
    if (!file) throw std::runtime_error("cannot write " + temporary.string());
    const json header = {{"format", "ellidft-checkpoint-1"},
                         {"fingerprint", checkpoint.fingerprint},
                         {"scf_converged", checkpoint.scf_converged},
                         {"n", checkpoint.density.rows()},
                         {"response", checkpoint.response}};
    file << header.dump() << "\n";
    file.write(reinterpret_cast<const char*>(checkpoint.density.data()),
               static_cast<std::streamsize>(sizeof(double)) * checkpoint.density.size());
    if (!file) throw std::runtime_error("cannot write " + temporary.string());
  }
  std::filesystem::rename(temporary, path);
}

// False if there is no usable file; `reason` then says why.
bool read_checkpoint(const std::filesystem::path& path, const std::string& fingerprint, Eigen::Index n,
                     Checkpoint& checkpoint, std::string& reason) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    reason = "no checkpoint file";
    return false;
  }
  try {
    std::string line;
    std::getline(file, line);
    const json header = json::parse(line);
    if (header.at("format").get<std::string>() != "ellidft-checkpoint-1") {
      reason = "unknown checkpoint format";
      return false;
    }
    if (header.at("fingerprint").get<std::string>() != fingerprint || header.at("n").get<Eigen::Index>() != n) {
      reason = "the checkpoint belongs to another molecule, basis or method";
      return false;
    }
    checkpoint.fingerprint = fingerprint;
    checkpoint.scf_converged = header.at("scf_converged").get<bool>();
    checkpoint.response = header.at("response");
    checkpoint.density.resize(n, n);
    file.read(reinterpret_cast<char*>(checkpoint.density.data()), static_cast<std::streamsize>(sizeof(double)) * n * n);
    if (!file || !checkpoint.density.allFinite()) {
      reason = "the checkpoint file is incomplete";
      return false;
    }
    return true;
  } catch (const std::exception& e) {
    reason = std::string("the checkpoint file cannot be read: ") + e.what();
    return false;
  }
}

Eigen::Matrix3d matrix_from_json(const json& rows) {
  Eigen::Matrix3d m;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) m(a, b) = rows.at(static_cast<std::size_t>(a)).at(static_cast<std::size_t>(b)).get<double>();
  return m;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: ellidft <control.json>\n";
      return EXIT_FAILURE;
    }
    const std::filesystem::path control_path = std::filesystem::absolute(argv[1]);
    std::ifstream in(control_path);
    if (!in) throw std::runtime_error("cannot read " + control_path.string());
    const json control = json::parse(in, nullptr, true, true);   // comments allowed
    check_keys(control, "control file", {"name", "geometry", "bohr", "charge", "basis", "aux_basis", "functional", "grid",
                                         "scf", "response", "threads", "output_directory", "restart"});
    const json scf_node = control.value("scf", json::object());
    const json response_node = control.value("response", json::object());
    check_keys(scf_node, "scf", {"energy_tolerance", "error_tolerance", "max_iterations", "guess", "overlap_threshold"});
    check_keys(response_node, "response", {"wavelengths_nm", "frequencies_hartree", "gamma_hartree", "states", "tda"});
    if (!control.contains("geometry") || !control.contains("basis"))
      throw std::runtime_error("control file: \"geometry\" and \"basis\" are required");

    const std::filesystem::path directory = control_path.parent_path();
    const std::filesystem::path output_directory = directory / control.value("output_directory", std::string("."));
    std::filesystem::create_directories(output_directory);
    const std::string name = control.value("name", control_path.stem().string());
    Protocol log(output_directory / (name + ".log"));
    const auto start = std::chrono::steady_clock::now();

    const int threads = control.value("threads", 0);
    if (threads > 0) set_num_threads(threads);

    log.line(format("elliDFT %s (commit %s, %s, %s)", ELLIDFT_VERSION, ELLIDFT_GIT_COMMIT, ELLIDFT_BUILD_TYPE, ELLIDFT_COMPILER));
    log.line(format("libraries: Eigen %d.%d.%d, libxc %s, nlohmann-json %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION,
                    EIGEN_MINOR_VERSION, libxc_version(), NLOHMANN_JSON_VERSION_MAJOR, NLOHMANN_JSON_VERSION_MINOR,
                    NLOHMANN_JSON_VERSION_PATCH));
    log.line("started " + now() + ", " + std::to_string(max_threads()) + " thread(s)");
    log.line("control file " + control_path.string());
    log.line();

    // --- input -------------------------------------------------------------
    const std::string geometry_file = (directory / control.at("geometry").get<std::string>()).string();
    Molecule molecule = read_xyz_file(geometry_file, control.value("bohr", false) ? LengthUnit::Bohr : LengthUnit::Angstrom);
    molecule.charge = control.value("charge", 0);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(control.at("basis").get<std::string>()));
    Method method;
    method.functional = Functional::by_name(control.value("functional", std::string("CAM-B3LYP")));
    if (control.contains("grid")) {
      if (!control.at("grid").is_array() || control.at("grid").size() != 2)
        throw std::runtime_error("control file: \"grid\" expects [radial, angular]");
      method.grid.n_radial = control.at("grid").at(0).get<int>();
      method.grid.n_angular = control.at("grid").at(1).get<int>();
    }
    const std::string aux_name = control.value("aux_basis", std::string());

    ScfOptions scf_options;
    scf_options.energy_tolerance = scf_node.value("energy_tolerance", scf_options.energy_tolerance);
    scf_options.error_tolerance = scf_node.value("error_tolerance", scf_options.error_tolerance);
    scf_options.max_iterations = scf_node.value("max_iterations", scf_options.max_iterations);
    scf_options.overlap_threshold = scf_node.value("overlap_threshold", scf_options.overlap_threshold);
    const std::string guess = scf_node.value("guess", std::string("atomic"));
    if (guess == "atomic") scf_options.initial_guess = InitialGuess::AtomicDensities;
    else if (guess == "core") scf_options.initial_guess = InitialGuess::CoreHamiltonian;
    else throw std::runtime_error("scf.guess expects \"atomic\" or \"core\"");

    std::vector<double> frequencies;
    if (response_node.contains("frequencies_hartree"))
      for (double w : response_node.at("frequencies_hartree").get<std::vector<double>>()) frequencies.push_back(w);
    if (response_node.contains("wavelengths_nm"))
      for (double nm : response_node.at("wavelengths_nm").get<std::vector<double>>()) {
        if (!(nm > 0.0)) throw std::runtime_error("response.wavelengths_nm: wavelengths must be positive");
        frequencies.push_back(optics::ev_nm / nm / units::hartree_in_ev);   // specification, section 8
      }
    const double gamma = response_node.value("gamma_hartree", 0.0);
    const int n_states = response_node.value("states", 0);
    const bool tda = response_node.value("tda", false);

    log.line("molecule   " + molecule.name);
    log.line(format("           %d atoms, charge %d, geometry %s", static_cast<int>(molecule.n_atoms()), molecule.charge, geometry_file.c_str()));
    log.line(format("basis      %s: %d functions, %d shells%s", basis.name().c_str(), basis.n_functions(),
                    static_cast<int>(basis.shells().size()), aux_name.empty() ? "" : (", density fitting with " + aux_name).c_str()));
    log.line(format("method     %s%s", method.functional.name.c_str(),
                    method.functional.has_density_part()
                        ? format(", grid (%d, %d)", method.grid.n_radial, method.grid.n_angular).c_str() : ""));
    log.line(format("scf        energy %.1e Hartree, commutator %.1e, at most %d iterations, start \"%s\"",
                    scf_options.energy_tolerance, scf_options.error_tolerance, scf_options.max_iterations, guess.c_str()));
    log.line(format("response   %d frequencies, damping %.3e Hartree, %d states%s", static_cast<int>(frequencies.size()), gamma,
                    n_states, tda ? " (Tamm-Dancoff)" : ""));
    log.line();

    // --- checkpoint ----------------------------------------------------------
    // What the stored data belong to: geometry, basis, method, thresholds that change the result.
    std::ostringstream identity;
    identity << std::setprecision(12) << "charge " << molecule.charge << " basis " << basis.name() << " "
             << basis.n_functions() << " aux " << aux_name << " functional " << method.functional.name << " grid "
             << method.grid.n_radial << "," << method.grid.n_angular << " overlap " << scf_options.overlap_threshold;
    for (const Atom& a : molecule.atoms) identity << " " << a.Z << " " << a.r[0] << " " << a.r[1] << " " << a.r[2];
    const std::filesystem::path checkpoint_path = output_directory / (name + ".chk");
    Checkpoint checkpoint;
    checkpoint.fingerprint = identity.str();
    bool restarted = false;
    if (control.value("restart", true)) {
      Checkpoint stored;
      std::string reason;
      if (read_checkpoint(checkpoint_path, checkpoint.fingerprint, basis.n_functions(), stored, reason)) {
        checkpoint = stored;
        restarted = true;
        scf_options.start_density = checkpoint.density;
        log.line(format("restart    from %s: density of %s SCF, %d finished response result(s)", checkpoint_path.c_str(),
                        checkpoint.scf_converged ? "a converged" : "an unfinished", static_cast<int>(checkpoint.response.size())));
      } else {
        log.line("restart    not used: " + reason);
        checkpoint.response = json::object();
      }
      log.line();
    }
    checkpoint.scf_converged = false;
    scf_options.on_iteration = [&](const ScfIteration&, const Eigen::MatrixXd& density) {
      checkpoint.density = density;
      write_checkpoint(checkpoint_path, checkpoint);
    };

    // --- ground state --------------------------------------------------------
    std::unique_ptr<JKBuilder> jk;
    if (aux_name.empty()) {
      jk = std::make_unique<DirectJK>(basis);
    } else {
      const BasisSet aux = BasisSet::build(molecule, BasisLibrary::load(aux_name));
      jk = std::make_unique<DensityFittedJK>(basis, aux);
    }
    const ScfResult scf = run_scf(molecule, basis, method, *jk, scf_options);
    log.line("self-consistent field");
    log.line("  iteration   total energy / Hartree      change        commutator");
    for (const ScfIteration& it : scf.history)
      log.line(format("  %6d   %22.12f   %11.3e   %11.3e", it.iteration, it.energy, it.energy_change, it.error));
    if (!scf.converged) {
      log.line("  NOT CONVERGED");
      throw std::runtime_error("the SCF did not converge");
    }
    log.line(format("  converged after %d iterations; %d orbitals (%d occupied), smallest overlap eigenvalue %.2e",
                    scf.iterations, scf.n_orbitals, scf.n_occupied, scf.overlap_min_eigenvalue));
    checkpoint.density = scf.density;
    checkpoint.scf_converged = true;
    write_checkpoint(checkpoint_path, checkpoint);
    log.line(format("  total energy %.10f Hartree", scf.energy_total));
    log.line(format("  one-electron %.10f, two-electron %.10f, exchange-correlation %.10f, nuclear repulsion %.10f",
                    scf.energy_one_electron, scf.energy_two_electron, scf.energy_xc, scf.energy_nuclear_repulsion));
    if (method.functional.has_density_part())
      log.line(format("  electrons on the grid %.8f (error per electron %.1e)", scf.electrons_on_grid,
                      std::abs(scf.electrons_on_grid - 2.0 * scf.n_occupied) / (2.0 * scf.n_occupied)));
    log.line(format("  dipole moment %.6f %.6f %.6f e Bohr", scf.dipole[0], scf.dipole[1], scf.dipole[2]));
    log.line(format("  time %.2f s: one-electron %.2f, start density %.2f, Coulomb and exchange %.2f, exchange-correlation %.2f, "
                    "diagonalisation %.2f", scf.timings.total, scf.timings.one_electron, scf.timings.initial_guess,
                    scf.timings.coulomb_exchange, scf.timings.exchange_correlation, scf.timings.diagonalisation));
    log.line();

    // --- response ------------------------------------------------------------
    const auto response_start = std::chrono::steady_clock::now();
    const ResponseOperator op(molecule, basis, method, scf, *jk);
    ordered_json symbols = ordered_json::array(), geometry = ordered_json::array();
    for (const Atom& a : molecule.atoms) {
      symbols.push_back(element_symbol(a.Z));
      geometry.push_back({a.r[0], a.r[1], a.r[2]});
    }
    ordered_json out;
    out["molecule"] = {{"name", molecule.name}, {"charge", molecule.charge}, {"symbols", symbols}, {"geometry_bohr", geometry}};
    out["method"] = {{"functional", method.functional.name},
                     {"basis", basis.name()},
                     {"aux_basis", aux_name},
                     {"grid", method.functional.has_density_part()
                                  ? std::to_string(method.grid.n_radial) + "," + std::to_string(method.grid.n_angular)
                                  : std::string()},
                     {"response", tda ? "TDA" : "RPA"}};
    out["scf"] = {{"energy_hartree", scf.energy_total}, {"dipole_au", {scf.dipole[0], scf.dipole[1], scf.dipole[2]}}};

    log.line("linear response");
    // Polarisability at one frequency: from the checkpoint if it is there, otherwise computed and stored.
    auto polarizability_with_checkpoint = [&](double omega, double damping, bool& from_checkpoint) {
      const std::string key = format("%.17g/%.17g", omega, damping);
      Polarizability alpha;
      from_checkpoint = restarted && checkpoint.response.contains(key);
      if (from_checkpoint) {
        const json& entry = checkpoint.response.at(key);
        alpha.omega = omega;
        alpha.gamma = damping;
        alpha.real = matrix_from_json(entry.at("re"));
        alpha.imag = matrix_from_json(entry.at("im"));
        alpha.iterations = entry.at("iterations").get<int>();
        alpha.converged = true;
        return alpha;
      }
      alpha = polarizability(op, omega, damping);
      if (alpha.converged) {
        checkpoint.response[key] = {{"re", matrix_json(alpha.real)}, {"im", matrix_json(alpha.imag)}, {"iterations", alpha.iterations}};
        write_checkpoint(checkpoint_path, checkpoint);
      }
      return alpha;
    };
    bool from_checkpoint = false;
    const Polarizability alpha0 = polarizability_with_checkpoint(0.0, 0.0, from_checkpoint);
    if (!alpha0.converged) throw std::runtime_error("the static response equations did not converge");
    out["polarizability_static_au"] = matrix_json(alpha0.real);
    log.line(format("  static polarisability, diagonal %.6f %.6f %.6f, mean %.6f a.u. (%d iterations%s)", alpha0.real(0, 0),
                    alpha0.real(1, 1), alpha0.real(2, 2), alpha0.real.trace() / 3.0, alpha0.iterations,
                    from_checkpoint ? ", from the checkpoint" : ""));
    out["polarizability_dynamic"] = ordered_json::array();
    for (double omega : frequencies) {
      const Polarizability alpha = polarizability_with_checkpoint(omega, gamma, from_checkpoint);
      if (!alpha.converged)
        throw std::runtime_error("the response equations did not converge at omega = " + std::to_string(omega));
      out["polarizability_dynamic"].push_back({{"omega_hartree", omega},
                                               {"gamma_hartree", gamma},
                                               {"alpha_re_au", matrix_json(alpha.real)},
                                               {"alpha_im_au", matrix_json(alpha.imag)}});
      log.line(format("  omega %.6f Hartree (%.1f nm): mean %.6f + %.6f i a.u. (%d iterations%s)", omega,
                      omega > 0.0 ? optics::ev_nm / (omega * units::hartree_in_ev) : 0.0, alpha.real.trace() / 3.0,
                      alpha.imag.trace() / 3.0, alpha.iterations, from_checkpoint ? ", from the checkpoint" : ""));
    }
    out["excitations"] = ordered_json::array();
    if (n_states > 0) {
      ExcitationOptions options;
      options.n_states = n_states;
      options.tamm_dancoff = tda;
      const ExcitationResult result = excitations(op, options);
      if (!result.converged) throw std::runtime_error("the excitation energies did not converge");
      log.line(format("  excitations (%d iterations):   energy / eV   wavelength / nm   oscillator strength", result.iterations));
      for (const Excitation& s : result.states) {
        out["excitations"].push_back({{"energy_hartree", s.energy},
                                      {"oscillator_strength", s.oscillator_strength},
                                      {"transition_dipole_au", {s.transition_dipole[0], s.transition_dipole[1], s.transition_dipole[2]}}});
        log.line(format("    %12.5f   %12.2f   %12.6f", s.energy * units::hartree_in_ev,
                        optics::ev_nm / (s.energy * units::hartree_in_ev), s.oscillator_strength));
      }
    }
    const double response_seconds = seconds_since(response_start);
    log.line(format("  time %.2f s, %d applications of the response operator", response_seconds, op.applications()));
    const auto* df = dynamic_cast<const DensityFittedJK*>(jk.get());
    if (df != nullptr)
      log.line(format("  density fitting: %d auxiliary functions, %.3f GB held, %.2f s to build", static_cast<int>(df->n_auxiliary()),
                      df->memory_bytes() / 1073741824.0, df->build_seconds()));
    log.line();

    out["run"] = {{"program", "ellidft"},
                  {"version", ELLIDFT_VERSION},
                  {"commit", ELLIDFT_GIT_COMMIT},
                  {"threads", max_threads()},
                  {"scf_iterations", scf.iterations},
                  {"scf_seconds", scf.timings.total},
                  {"response_seconds", response_seconds},
                  {"response_operator_applications", op.applications()},
                  {"total_seconds", seconds_since(start)}};
    if (df != nullptr)
      out["run"]["density_fitting"] = {{"auxiliary_functions", df->n_auxiliary()},
                                       {"memory_gb", df->memory_bytes() / 1073741824.0},
                                       {"build_seconds", df->build_seconds()}};
    const std::filesystem::path result_path = output_directory / (name + ".json");
    std::ofstream result_file(result_path);
    if (!result_file) throw std::runtime_error("cannot write " + result_path.string());
    result_file << out.dump(2) << "\n";
    log.line("results written to " + result_path.string());
    log.line(format("finished %s, total time %.2f s", now().c_str(), seconds_since(start)));
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
