// ellidft-nonlinear: non-linear optical quantities of a closed-shell molecule
// from real-time propagation, for disordered materials (photopolymers):
// optical Kerr effect (orientational average, optionally n2 of a film),
// two-photon absorption (strengths per state, cross-section spectrum) and
// second-harmonic generation along a direction.
//
// Usage: ellidft-nonlinear <control.json>      (docs/Nichtlineare-Optik.md)
//
//   {
//     "geometry": "molecule.xyz", "basis": "def2-SVPD", "functional": "CAM-B3LYP",
//     "grid": [75, 302], "aux_basis": "def2-universal-JKFIT", "charge": 0,
//     "time_step": 0.2, "fock_tolerance": 1e-11,
//     "kerr": {"wavelengths_nm": [780], "amplitude": 0.01, "ramp_time": 300, "periods": 2,
//              "film": {"density_g_cm3": 1.18, "molar_mass_g_mol": 100.12, "refractive_index": 1.49}},
//     "two_photon": {"states": 6, "amplitude": 0.01, "duration": 300, "half_width_ev": 0.1,
//                    "wavelengths_nm": {"from": 600, "to": 1000, "step": 10}},
//     "second_harmonic": {"wavelengths_nm": [1064], "direction": [0, 0, 1],
//                         "amplitude": 0.01, "ramp_time": 300, "periods": 2}
//   }
//
// Every task is optional. Results as JSON on standard output, a protocol on
// standard error. File names are relative to the control file.
//
// Finished parts (every wavelength of "kerr" and "second_harmonic", the
// strengths of "two_photon") are written to <control file>.partial.json at
// once. A later run with the same molecule, method and settings takes them
// from there instead of computing them again; "restart": false switches
// this off.

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/parallel.hpp"
#include "core/units.hpp"
#include "core/xyz.hpp"
#include "df/df_jk.hpp"
#include "optics/optics.hpp"
#include "realtime/nonlinear.hpp"
#include "response/excitations.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using nlohmann::json;
using nlohmann::ordered_json;

namespace {

void check_keys(const json& node, const std::string& where, std::initializer_list<const char*> allowed) {
  for (const auto& [key, unused] : node.items()) {
    bool known = key == "comment";
    for (const char* a : allowed) known = known || key == a;
    if (!known) throw std::runtime_error(where + ": unknown key \"" + key + "\"");
  }
}

// Photon energy in Hartree for a wavelength in nm (specification, section 8).
double photon_energy(double wavelength_nm) { return optics::ev_nm / wavelength_nm / units::hartree_in_ev; }

std::vector<double> wavelengths(const json& node, const std::string& where) {
  std::vector<double> values;
  if (node.is_array()) {
    values = node.get<std::vector<double>>();
  } else if (node.is_number()) {
    values.push_back(node.get<double>());
  } else {
    const double from = node.at("from").get<double>(), to = node.at("to").get<double>(), step = node.at("step").get<double>();
    if (!(step > 0.0) || to < from) throw std::runtime_error(where + ": need from <= to and step > 0");
    for (int i = 0; from + i * step <= to + 1e-9 * step; ++i) values.push_back(from + i * step);
  }
  for (double v : values)
    if (!(v > 0.0)) throw std::runtime_error(where + ": wavelengths must be positive");
  if (values.empty()) throw std::runtime_error(where + ": no wavelengths");
  return values;
}

std::string format_energy(double energy) {
  char buffer[40];
  std::snprintf(buffer, sizeof buffer, "%.10f ", energy);
  return buffer;
}

HarmonicFieldOptions harmonic_options(const json& node, const RealTimeOptions& propagation) {
  HarmonicFieldOptions options;
  options.amplitude = node.value("amplitude", 0.01);
  options.ramp_time = node.value("ramp_time", 300.0);
  options.n_periods = node.value("periods", 2);
  options.propagation = propagation;
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: ellidft-nonlinear <control.json>\n";
      return EXIT_FAILURE;
    }
    std::ifstream in(argv[1]);
    if (!in) throw std::runtime_error(std::string("cannot read ") + argv[1]);
    const json control = json::parse(in, nullptr, true, true);
    check_keys(control, "control file", {"geometry", "basis", "functional", "grid", "aux_basis", "charge", "bohr",
                                         "time_step", "fock_tolerance", "kerr", "two_photon", "second_harmonic", "restart"});
    const std::filesystem::path directory = std::filesystem::absolute(argv[1]).parent_path();

    Molecule molecule = read_xyz_file((directory / control.at("geometry").get<std::string>()).string(),
                                      control.value("bohr", false) ? LengthUnit::Bohr : LengthUnit::Angstrom);
    molecule.charge = control.value("charge", 0);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(control.at("basis").get<std::string>()));
    Method method;
    method.functional = Functional::by_name(control.value("functional", std::string("CAM-B3LYP")));
    if (control.contains("grid")) {
      method.grid.n_radial = control.at("grid").at(0).get<int>();
      method.grid.n_angular = control.at("grid").at(1).get<int>();
    }
    const std::string aux_name = control.value("aux_basis", std::string());
    std::unique_ptr<JKBuilder> jk;
    if (aux_name.empty()) {
      jk = std::make_unique<DirectJK>(basis);
    } else {
      const BasisSet aux = BasisSet::build(molecule, BasisLibrary::load(aux_name));
      jk = std::make_unique<DensityFittedJK>(basis, aux);
    }
    ScfOptions scf_options;
    scf_options.energy_tolerance = 1e-12;
    scf_options.error_tolerance = 1e-10;
    const ScfResult scf = run_scf(molecule, basis, method, *jk, scf_options);
    if (!scf.converged) throw std::runtime_error("the SCF did not converge");
    const RealTimePropagator propagator(molecule, basis, method, scf, *jk);

    RealTimeOptions propagation;
    propagation.time_step = control.value("time_step", 0.2);
    propagation.fock_tolerance = control.value("fock_tolerance", 1e-11);
    propagation.max_fock_builds_per_step = 20;

    ordered_json out;
    out["program"] = "ellidft-nonlinear";
    out["molecule"] = {{"name", molecule.name}, {"charge", molecule.charge}, {"atoms", molecule.n_atoms()}};
    out["method"] = {{"functional", method.functional.name}, {"basis", basis.name()}, {"aux_basis", aux_name},
                     {"time_step_au", propagation.time_step}};
    out["scf"] = {{"energy_hartree", scf.energy_total}};
    out["conventions"] = {
        {"expansion", "mu = mu0 + alpha F + 1/2 beta F F + 1/6 gamma F F F; frequency-dependent quantities reduce to these for omega -> 0"},
        {"averages", "over all orientations of the molecule, light linearly polarised, all fields parallel"},
        {"two_photon", "cross-section per absorbed photon, one beam, Lorentzian line; sigma = 4 pi^3 alpha^2 omega^2 strength g(2 omega)"},
        {"units", "atomic units unless the name says otherwise"}};
    std::cerr << std::setprecision(6);

    // --- finished parts of an earlier run -------------------------------------
    // Key of a part: what it was computed for (geometry, method, step) and the settings of its task.
    std::ostringstream identity;
    identity << std::setprecision(12) << "charge " << molecule.charge << " basis " << basis.name() << " aux " << aux_name
             << " functional " << method.functional.name << " grid " << method.grid.n_radial << "," << method.grid.n_angular
             << " step " << propagation.time_step << " tolerance " << propagation.fock_tolerance;
    for (const Atom& a : molecule.atoms) identity << " " << a.Z << " " << a.r[0] << " " << a.r[1] << " " << a.r[2];
    const std::string common_key = identity.str();
    std::filesystem::path partial_path = std::filesystem::absolute(argv[1]);
    partial_path += ".partial.json";
    json partial = json::object();
    if (control.value("restart", true)) {
      std::ifstream stored(partial_path);
      if (stored) {
        try {
          partial = json::parse(stored);
          if (!partial.is_object()) partial = json::object();
        } catch (const std::exception&) {
          partial = json::object();   // unreadable: start again
        }
      }
    }
    auto task_key = [&](const std::string& task, json settings, const std::string& item) {
      settings.erase("wavelengths_nm");   // the wavelengths are items, not settings
      settings.erase("comment");
      return task + " | " + common_key + " | " + settings.dump() + " | " + item;
    };
    auto store = [&](const std::string& key, const json& value) {
      partial[key] = value;
      std::filesystem::path temporary = partial_path;
      temporary += ".tmp";
      {
        std::ofstream file(temporary);
        if (!file) throw std::runtime_error("cannot write " + temporary.string());
        file << partial.dump() << "\n";
      }
      std::filesystem::rename(temporary, partial_path);
    };

    if (control.contains("kerr")) {
      const json& node = control.at("kerr");
      check_keys(node, "kerr", {"wavelengths_nm", "amplitude", "ramp_time", "periods", "film"});
      ordered_json rows = ordered_json::array();
      for (double nm : wavelengths(node.at("wavelengths_nm"), "kerr.wavelengths_nm")) {
        HarmonicFieldOptions options = harmonic_options(node, propagation);
        options.omega = photon_energy(nm);
        json film_free = node;
        film_free.erase("film");   // the film data do not enter the propagation
        const std::string key = task_key("kerr", film_free, std::to_string(nm));
        ordered_json row;
        const bool finished = partial.contains(key);
        if (finished) {
          row = partial.at(key);
        } else {
          const KerrAverage kerr = kerr_average(propagator, options);
          row = {{"wavelength_nm", nm},
                 {"omega_hartree", options.omega},
                 {"gamma_average_au", kerr.gamma},
                 {"gamma_quadrature_au", kerr.gamma_quadrature},
                 {"gamma_by_direction_au", kerr.by_direction},
                 {"largest_fit_residual", kerr.largest_residual},
                 {"steps_converged", kerr.steps_converged},
                 {"fock_builds", kerr.fock_builds}};
          store(key, row);
        }
        struct {
          double gamma, largest_residual;
        } kerr{row.at("gamma_average_au").get<double>(), row.at("largest_fit_residual").get<double>()};
        std::cerr << "Kerr effect at " << nm << " nm: gamma(-w; w, w, -w) averaged = " << kerr.gamma << " a.u., fit residual "
                  << kerr.largest_residual << (finished ? "  (from the earlier run)" : "") << "\n";
        if (kerr.largest_residual > 0.05)
          std::cerr << "  WARNING: the harmonic fit leaves more than 5 % (resonance near, or switch-on too fast)\n";
        if (node.contains("film")) {
          const json& film = node.at("film");
          check_keys(film, "kerr.film", {"density_g_cm3", "molar_mass_g_mol", "refractive_index"});
          const double N = optics::number_density(film.at("density_g_cm3").get<double>(), film.at("molar_mass_g_mol").get<double>());
          const NonlinearIndex index = nonlinear_index(kerr.gamma, N, film.at("refractive_index").get<double>());
          row["film"] = {{"number_density_per_A3", N},
                         {"local_field_factor", index.local_field_factor},
                         {"n2_m2_per_W", index.n2_m2_per_W},
                         {"n2_cm2_per_W", index.n2_cm2_per_W},
                         {"note", "electronic part, Lorentz local field; the molecule of this calculation is taken as the unit"}};
          std::cerr << "  film: n2 = " << index.n2_cm2_per_W << " cm^2/W\n";
        }
        rows.push_back(row);
      }
      out["kerr"] = rows;
    }

    if (control.contains("two_photon")) {
      const json& node = control.at("two_photon");
      check_keys(node, "two_photon", {"states", "amplitude", "duration", "half_width_ev", "wavelengths_nm", "margin_hartree"});
      const ResponseOperator op(molecule, basis, method, scf, *jk);
      ExcitationOptions excitation_options;
      excitation_options.n_states = node.value("states", 6);
      const ExcitationResult excited = excitations(op, excitation_options);
      if (!excited.converged) throw std::runtime_error("the excitation energies did not converge");
      double first_allowed = 0.0;
      for (const Excitation& s : excited.states)
        if (first_allowed == 0.0 && s.oscillator_strength > 1e-4) first_allowed = s.energy;
      if (first_allowed == 0.0) first_allowed = excited.states.front().energy;
      // The method needs photons below every one-photon line: states with half
      // their energy closer than the margin to the first allowed excitation are left out.
      const double margin = node.value("margin_hartree", 0.05);
      TwoPhotonStrengthOptions options;
      options.amplitude = node.value("amplitude", 0.01);
      options.duration = node.value("duration", 300.0);
      options.propagation = propagation;
      options.propagation.propagator = Propagator::Composition4;
      ordered_json skipped = ordered_json::array();
      for (const Excitation& s : excited.states) {
        if (0.5 * s.energy < first_allowed - margin) options.state_energies.push_back(s.energy);
        else skipped.push_back(s.energy);
      }
      ordered_json result = {{"first_allowed_excitation_hartree", first_allowed},
                             {"states_left_out_hartree", skipped},
                             {"half_width_ev", node.value("half_width_ev", 0.1)}};
      if (options.state_energies.empty()) {
        std::cerr << "two-photon absorption: no state with half its energy below the first allowed excitation\n";
        result["states"] = ordered_json::array();
      } else {
        std::string energies;
        for (double e : options.state_energies) energies += format_energy(e);
        const std::string key = task_key("two_photon", node, energies);
        TwoPhotonStrengths strengths;
        if (partial.contains(key)) {
          const json& entry = partial.at(key);
          for (const json& s : entry.at("states"))
            strengths.states.push_back({s.at(0).get<double>(), s.at(1).get<double>()});
          strengths.steps_converged = entry.at("steps_converged").get<bool>();
          strengths.fock_builds = entry.at("fock_builds").get<int>();
          std::cerr << "two-photon absorption: strengths from the earlier run\n";
        } else {
          strengths = two_photon_strengths(propagator, options);
          json list = json::array();
          for (const TwoPhotonState& s : strengths.states) list.push_back({s.energy, s.strength});
          store(key, {{"states", list}, {"steps_converged", strengths.steps_converged}, {"fock_builds", strengths.fock_builds}});
        }
        ordered_json states = ordered_json::array();
        for (const TwoPhotonState& s : strengths.states) {
          states.push_back({{"energy_hartree", s.energy},
                            {"energy_ev", s.energy * units::hartree_in_ev},
                            {"two_photon_wavelength_nm", optics::ev_nm / (0.5 * s.energy * units::hartree_in_ev)},
                            {"strength_average_au", s.strength}});
          std::cerr << "two-photon absorption: state at " << s.energy * units::hartree_in_ev << " eV, strength " << s.strength
                    << " a.u." << (s.strength < 0.0 ? "  WARNING: negative, not a valid result" : "") << "\n";
        }
        result["states"] = states;
        result["steps_converged"] = strengths.steps_converged;
        result["fock_builds"] = strengths.fock_builds;
        if (node.contains("wavelengths_nm")) {
          const double half_width = node.value("half_width_ev", 0.1) / units::hartree_in_ev;
          ordered_json spectrum = ordered_json::array();
          for (double nm : wavelengths(node.at("wavelengths_nm"), "two_photon.wavelengths_nm"))
            spectrum.push_back({nm, two_photon_cross_section_gm(strengths.states, photon_energy(nm), half_width)});
          result["cross_section_gm"] = spectrum;
          result["cross_section_columns"] = {"wavelength_nm", "sigma2_GM"};
        }
      }
      out["two_photon"] = result;
    }

    if (control.contains("second_harmonic")) {
      const json& node = control.at("second_harmonic");
      check_keys(node, "second_harmonic", {"wavelengths_nm", "direction", "amplitude", "ramp_time", "periods"});
      const std::vector<double> d = node.at("direction").get<std::vector<double>>();
      if (d.size() != 3) throw std::runtime_error("second_harmonic.direction: expected three numbers");
      const double norm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      if (!(norm > 0.0)) throw std::runtime_error("second_harmonic.direction: zero vector");
      const Vec3 u{d[0] / norm, d[1] / norm, d[2] / norm};
      ordered_json rows = ordered_json::array();
      for (double nm : wavelengths(node.at("wavelengths_nm"), "second_harmonic.wavelengths_nm")) {
        HarmonicFieldOptions options = harmonic_options(node, propagation);
        options.omega = photon_energy(nm);
        options.optical = u;
        const std::string key = task_key("second_harmonic", node, std::to_string(nm));
        ordered_json row;
        const bool finished = partial.contains(key);
        if (finished) {
          row = partial.at(key);
        } else {
          const HarmonicResponse r = harmonic_response(propagator, options);
          const Vec3& c = r.cosine[1][2];
          row = {{"wavelength_nm", nm},
                 {"omega_hartree", options.omega},
                 {"beta_vector_au", {4.0 * c[0], 4.0 * c[1], 4.0 * c[2]}},
                 {"beta_along_direction_au", 4.0 * (c[0] * u[0] + c[1] * u[1] + c[2] * u[2])},
                 {"fit_residual", r.residual[1]},
                 {"steps_converged", r.steps_converged}};
          store(key, row);
        }
        rows.push_back(row);
        std::cerr << "second harmonic at " << nm << " nm: beta(-2w; w, w) along the direction = "
                  << row.at("beta_along_direction_au").get<double>() << " a.u." << (finished ? "  (from the earlier run)" : "")
                  << "\n";
      }
      out["second_harmonic"] = {{"direction", {u[0], u[1], u[2]}},
                                {"note", "molecular quantity for fields along the direction; in a disordered material "
                                         "second-harmonic generation vanishes in the bulk"},
                                {"results", rows}};
    }

    out["run"] = {{"threads", max_threads()}, {"scf_iterations", scf.iterations}};
    std::cout << out.dump(2) << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
