// ellidft-kerr-film: optical Kerr effect of a polymer film from an oligomer
// series. The second hyperpolarisability per repeat unit is the difference of
// consecutive oligomers, as for the polarisability (specification, section 8,
// "Vom Oligomer zu n und k"); n2 follows with the Lorentz local field
// (nonlinear_index in src/realtime/nonlinear.hpp).
//
// Usage: ellidft-kerr-film --density <g/cm^3> --molar-mass <g/mol>
//          (--refractive-index <n0> | --optics <output of ellidft-optics>)
//          [--name <text>] <1 unit>.json <2 units>.json [<3 units>.json ...]
//
// The input files are outputs of ellidft-nonlinear with the task "kerr" for
// the oligomers with 1, 2, 3, ... repeat units, in this order, at the same
// wavelengths. With --optics the linear index is taken from that file at the
// same wavelengths (they must be present there); otherwise it is the given
// constant. Prints JSON; n2 uses the increment of the two longest oligomers.

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "optics/optics.hpp"
#include "realtime/nonlinear.hpp"

using nlohmann::json;
using nlohmann::ordered_json;

namespace {

json read_json(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  try {
    return json::parse(in);
  } catch (const std::exception& e) {
    throw std::runtime_error(path + ": " + e.what());
  }
}

struct Oligomer {
  std::string file;
  std::vector<double> wavelengths, gamma, residual;
};

Oligomer read_oligomer(const std::string& file) {
  const json data = read_json(file);
  if (!data.contains("kerr")) throw std::runtime_error(file + ": no \"kerr\" results (output of ellidft-nonlinear expected)");
  Oligomer o;
  o.file = file;
  for (const json& row : data.at("kerr")) {
    if (!row.value("steps_converged", true)) throw std::runtime_error(file + ": a propagation did not converge");
    o.wavelengths.push_back(row.at("wavelength_nm").get<double>());
    o.gamma.push_back(row.at("gamma_average_au").get<double>());
    o.residual.push_back(row.value("largest_fit_residual", 0.0));
  }
  if (o.wavelengths.empty()) throw std::runtime_error(file + ": no wavelengths");
  return o;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace ellidft;
  try {
    double density = 0.0, molar_mass = 0.0, n0 = 0.0;
    std::string name, optics_file;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      auto value = [&]() -> std::string {
        if (i + 1 >= argc) throw std::invalid_argument("option " + arg + " needs a value");
        return argv[++i];
      };
      if (arg == "--density") density = std::stod(value());
      else if (arg == "--molar-mass") molar_mass = std::stod(value());
      else if (arg == "--refractive-index") n0 = std::stod(value());
      else if (arg == "--optics") optics_file = value();
      else if (arg == "--name") name = value();
      else if (!arg.empty() && arg[0] == '-') throw std::invalid_argument("unknown option " + arg);
      else files.push_back(arg);
    }
    const bool have_index = n0 > 0.0;
    if (files.size() < 2 || !(density > 0.0) || !(molar_mass > 0.0) || have_index == !optics_file.empty()) {
      std::cerr << "usage: ellidft-kerr-film --density <g/cm^3> --molar-mass <g/mol> "
                   "(--refractive-index <n0> | --optics <file>) [--name <text>] <1 unit>.json <2 units>.json [...]\n";
      return EXIT_FAILURE;
    }

    std::vector<Oligomer> oligomers;
    for (const std::string& file : files) oligomers.push_back(read_oligomer(file));
    for (const Oligomer& o : oligomers)
      if (o.wavelengths != oligomers.front().wavelengths)
        throw std::runtime_error(o.file + ": wavelengths differ from " + oligomers.front().file);
    const std::size_t n_oligomers = oligomers.size(), n_wavelengths = oligomers.front().wavelengths.size();

    // linear index at every wavelength
    std::vector<double> index(n_wavelengths, n0);
    if (!optics_file.empty()) {
      const json optics_data = read_json(optics_file);
      for (std::size_t w = 0; w < n_wavelengths; ++w) {
        bool found = false;
        for (const json& row : optics_data.at("optical_constants"))
          if (std::abs(row.at("wavelength_nm").get<double>() - oligomers.front().wavelengths[w]) < 1e-6) {
            index[w] = row.at("n").get<double>();
            found = true;
          }
        if (!found)
          throw std::runtime_error(optics_file + ": no entry at " + std::to_string(oligomers.front().wavelengths[w]) + " nm");
      }
    }

    const double N = optics::number_density(density, molar_mass);
    ordered_json out;
    out["material"] = {{"name", name},
                       {"density_g_cm3", density},
                       {"molar_mass_repeat_unit_g_mol", molar_mass},
                       {"number_density_per_A3", N}};
    out["model"] = {{"increment", "gamma(m) - gamma(m-1), orientational average of gamma(-w; w, w, -w)"},
                    {"local_field", "Lorentz, (n0^2 + 2) / 3 for every field"},
                    {"scope", "electronic part only; no thermal contribution, no reorientation"},
                    {"linear_index", optics_file.empty() ? "constant given on the command line" : "from " + optics_file}};
    ordered_json rows = ordered_json::array();
    for (std::size_t w = 0; w < n_wavelengths; ++w) {
      ordered_json series = ordered_json::array();
      std::vector<double> increments;
      double largest_residual = 0.0;
      for (std::size_t m = 0; m < n_oligomers; ++m) largest_residual = std::max(largest_residual, oligomers[m].residual[w]);
      for (std::size_t m = 1; m < n_oligomers; ++m) {
        increments.push_back(oligomers[m].gamma[w] - oligomers[m - 1].gamma[w]);
        series.push_back({{"units", m + 1}, {"gamma_increment_au", increments.back()}});
      }
      const NonlinearIndex film = nonlinear_index(increments.back(), N, index[w]);
      ordered_json row = {{"wavelength_nm", oligomers.front().wavelengths[w]},
                          {"gamma_increment_series", series},
                          {"gamma_per_unit_au", increments.back()},
                          {"refractive_index", index[w]},
                          {"local_field_factor", film.local_field_factor},
                          {"n2_m2_per_W", film.n2_m2_per_W},
                          {"n2_cm2_per_W", film.n2_cm2_per_W},
                          {"largest_fit_residual", largest_residual}};
      if (increments.size() >= 2) {
        const double last = increments.back(), previous = increments[increments.size() - 2];
        row["gamma_increment_change_percent"] = 100.0 * std::abs(last - previous) / std::abs(last);
      }
      rows.push_back(row);
      std::cerr << "Kerr effect of the film at " << oligomers.front().wavelengths[w] << " nm: gamma per unit " << increments.back()
                << " a.u., n2 = " << film.n2_cm2_per_W << " cm^2/W\n";
    }
    out["kerr"] = rows;
    std::cout << out.dump(2) << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
