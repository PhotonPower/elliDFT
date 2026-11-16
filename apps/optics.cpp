// ellidft-optics: from the response of an oligomer series to n(lambda) and
// k(lambda) (specification, section 8).
//
// Usage: ellidft-optics --density <g/cm^3> --molar-mass <g/mol>
//          [--order-parameter <S>] [--frame lab|principal] [--name <text>] [--table <file>]
//          <1 unit>.json <2 units>.json [<3 units>.json ...]
//
// The input files are outputs of ellidft-response for the oligomers with
// 1, 2, 3, ... repeat units, in this order, all with the same frequencies
// and damping. The polarisability per repeat unit is the difference of
// consecutive oligomers; the optical constants use the increment of the two
// longest oligomers.
//
// Chain axis (--frame):
//   lab        (default) the z axis of the geometries is the chain axis. The
//              oligomers must share one orientation, as the files in
//              data/geometries do.
//   principal  every oligomer gets its own frame from the principal axes of
//              its atoms (z = largest extent). Only the mean of the increment
//              is meaningful then; parallel and perpendicular components are
//              not comparable between oligomers of different length.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "core/elements.hpp"
#include "ellipsometry/dispersion_fit.hpp"
#include "optics/optics.hpp"

namespace {

using Complex = std::complex<double>;
using nlohmann::json;
using nlohmann::ordered_json;

Eigen::Matrix3d read_matrix(const json& rows) {
  Eigen::Matrix3d m;
  for (int k = 0; k < 3; ++k)
    for (int l = 0; l < 3; ++l) m(k, l) = rows.at(static_cast<std::size_t>(k)).at(static_cast<std::size_t>(l)).get<double>();
  return m;
}

struct Oligomer {
  std::string file;
  ellidft::Molecule molecule;
  std::vector<double> omega, gamma;
  // Uniaxial components in the chain frame: static first, then one per frequency.
  ellidft::optics::UniaxialPolarizability static_alpha;
  std::vector<ellidft::optics::UniaxialPolarizability> dynamic_alpha;
};

Oligomer read_oligomer(const std::string& file, bool principal_frame) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("cannot open " + file);
  const json data = json::parse(in);
  Oligomer o;
  o.file = file;
  const json& symbols = data.at("molecule").at("symbols");
  const json& geometry = data.at("molecule").at("geometry_bohr");
  for (std::size_t i = 0; i < symbols.size(); ++i) {
    ellidft::Atom atom;
    atom.Z = ellidft::atomic_number(symbols[i].get<std::string>());
    for (std::size_t k = 0; k < 3; ++k) atom.r[k] = geometry[i].at(k).get<double>();
    o.molecule.atoms.push_back(atom);
  }
  const Eigen::Matrix3d frame =
      principal_frame ? ellidft::optics::chain_frame(o.molecule) : Eigen::Matrix3d(Eigen::Matrix3d::Identity());
  auto reduce = [&](const Eigen::Matrix3cd& alpha) {
    return ellidft::optics::uniaxial(ellidft::optics::rotate_into_frame(alpha, frame));
  };
  o.static_alpha = reduce(read_matrix(data.at("polarizability_static_au")).cast<Complex>());
  for (const json& entry : data.at("polarizability_dynamic")) {
    o.omega.push_back(entry.at("omega_hartree").get<double>());
    o.gamma.push_back(entry.at("gamma_hartree").get<double>());
    Eigen::Matrix3cd alpha = read_matrix(entry.at("alpha_re_au")).cast<Complex>();
    alpha += Complex(0.0, 1.0) * read_matrix(entry.at("alpha_im_au")).cast<Complex>();
    o.dynamic_alpha.push_back(reduce(alpha));
  }
  return o;
}

ordered_json complex_json(Complex value) { return {{"re", value.real()}, {"im", value.imag()}}; }

ellidft::optics::UniaxialPolarizability difference(const ellidft::optics::UniaxialPolarizability& a,
                                                   const ellidft::optics::UniaxialPolarizability& b) {
  return {a.parallel - b.parallel, a.perpendicular - b.perpendicular};
}

}  // namespace

int main(int argc, char** argv) {
  using namespace ellidft;
  try {
    double density = 0.0, molar_mass = 0.0, order_parameter = 0.0;
    std::string name, table_file, frame_name = "lab";
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      auto value = [&]() -> std::string {
        if (i + 1 >= argc) throw std::invalid_argument("option " + arg + " needs a value");
        return argv[++i];
      };
      if (arg == "--density") density = std::stod(value());
      else if (arg == "--molar-mass") molar_mass = std::stod(value());
      else if (arg == "--order-parameter") order_parameter = std::stod(value());
      else if (arg == "--name") name = value();
      else if (arg == "--frame") frame_name = value();
      else if (arg == "--table") table_file = value();
      else if (!arg.empty() && arg[0] == '-') throw std::invalid_argument("unknown option " + arg);
      else files.push_back(arg);
    }
    if (files.size() < 2 || !(density > 0.0) || !(molar_mass > 0.0)) {
      std::cerr << "usage: ellidft-optics --density <g/cm^3> --molar-mass <g/mol> [--order-parameter <S>] "
                   "[--frame lab|principal] [--name <text>] [--table <file>] <1 unit>.json <2 units>.json [...]\n";
      return EXIT_FAILURE;
    }

    if (frame_name != "lab" && frame_name != "principal")
      throw std::invalid_argument("--frame expects 'lab' or 'principal'");
    const bool principal_frame = frame_name == "principal";
    if (principal_frame && order_parameter != 0.0)
      throw std::invalid_argument("an order parameter needs a common chain axis: use --frame lab");
    std::vector<Oligomer> oligomers;
    for (const std::string& file : files) oligomers.push_back(read_oligomer(file, principal_frame));
    for (const Oligomer& o : oligomers)
      if (o.omega != oligomers.front().omega || o.gamma != oligomers.front().gamma)
        throw std::runtime_error(o.file + ": frequencies or damping differ from " + oligomers.front().file);

    const double N = optics::number_density(density, molar_mass);
    const std::size_t n_oligomers = oligomers.size();
    const std::size_t n_frequencies = oligomers.front().omega.size();

    ordered_json out;
    out["material"] = {{"name", name},
                       {"density_g_cm3", density},
                       {"molar_mass_repeat_unit_g_mol", molar_mass},
                       {"number_density_per_A3", N},
                       {"order_parameter", order_parameter}};
    out["model"] = {{"increment", "alpha(m) - alpha(m-1)"},
                    {"chain_axis", principal_frame ? "principal axis of each oligomer (components not comparable)"
                                                   : "z axis of the geometries, common to all oligomers"},
                    {"local_field", "Lorentz-Lorenz"},
                    {"anisotropy", order_parameter == 0.0
                                       ? "isotropic"
                                       : "APPROXIMATION: Lorentz-Lorenz applied separately to the ordinary and "
                                         "extraordinary component"}};

    // Static increments: the whole series, to show the convergence.
    ordered_json series = ordered_json::array();
    std::vector<double> static_means;
    for (std::size_t m = 1; m < n_oligomers; ++m) {
      const auto inc = difference(oligomers[m].static_alpha, oligomers[m - 1].static_alpha);
      static_means.push_back(inc.mean().real() * optics::angstrom3_per_au);
      series.push_back({{"units", m + 1},
                        {"mean_A3", inc.mean().real() * optics::angstrom3_per_au},
                        {"parallel_A3", inc.parallel.real() * optics::angstrom3_per_au},
                        {"perpendicular_A3", inc.perpendicular.real() * optics::angstrom3_per_au}});
    }
    out["static_increment_series"] = series;
    if (static_means.size() >= 2) {
      const double last = static_means.back(), previous = static_means[static_means.size() - 2];
      out["static_increment_change_percent"] = 100.0 * std::abs(last - previous) / std::abs(last);
    }
    out["static_increment_mean_A3"] = static_means.back();

    ordered_json rows = ordered_json::array();
    std::ofstream table;
    if (!table_file.empty()) {
      table.open(table_file);
      if (!table) throw std::runtime_error("cannot write " + table_file);
      table << "# " << (name.empty() ? "material" : name) << ": density " << density << " g/cm^3, repeat unit "
            << molar_mass << " g/mol, order parameter " << order_parameter << "\n";
      table << (order_parameter == 0.0 ? "# wavelength_nm n k\n"
                                       : "# wavelength_nm n k n_o k_o n_e k_e  (o/e: component-wise local field, approximation)\n");
    }
    for (std::size_t f = 0; f < n_frequencies; ++f) {
      const auto inc = difference(oligomers[n_oligomers - 1].dynamic_alpha[f], oligomers[n_oligomers - 2].dynamic_alpha[f]);
      Eigen::Matrix3cd unit = Eigen::Matrix3cd::Zero();
      unit(0, 0) = unit(1, 1) = inc.perpendicular;
      unit(2, 2) = inc.parallel;
      const double wavelength = optics::wavelength_nm(oligomers.front().omega[f]);
      const optics::OpticalConstants c = optics::optical_constants(wavelength, unit, N, order_parameter);
      ordered_json row = {{"wavelength_nm", wavelength},
                          {"omega_hartree", oligomers.front().omega[f]},
                          {"gamma_hartree", oligomers.front().gamma[f]},
                          {"increment_mean_A3", complex_json(inc.mean() * optics::angstrom3_per_au)},
                          {"n", c.isotropic.real()},
                          {"k", c.isotropic.imag()}};
      if (order_parameter != 0.0) {
        row["n_ordinary"] = c.ordinary.real();
        row["k_ordinary"] = c.ordinary.imag();
        row["n_extraordinary"] = c.extraordinary.real();
        row["k_extraordinary"] = c.extraordinary.imag();
      }
      rows.push_back(row);
      if (table.is_open()) {
        table << std::fixed << std::setprecision(3) << wavelength << " " << std::setprecision(6) << c.isotropic.real()
              << " " << c.isotropic.imag();
        if (order_parameter != 0.0)
          table << " " << c.ordinary.real() << " " << c.ordinary.imag() << " " << c.extraordinary.real() << " "
                << c.extraordinary.imag();
        table << "\n";
      }
    }
    out["optical_constants"] = rows;

    // Dispersion models through the computed index, for use as a material in ellidft-ellipsometry. Only
    // wavelengths without noticeable absorption take part.
    {
      std::vector<double> transparent_wavelengths, transparent_n;
      for (const auto& row : rows)
        if (row.at("k").get<double>() < 1e-4) {
          transparent_wavelengths.push_back(row.at("wavelength_nm").get<double>());
          transparent_n.push_back(row.at("n").get<double>());
        }
      ordered_json models = {{"points", transparent_wavelengths.size()},
                             {"note", "isotropic index, wavelengths with k < 1e-4; valid inside the computed range only"}};
      if (transparent_wavelengths.size() >= 3) {
        const auto cauchy = ellidft::ellipsometry::fit_cauchy(transparent_wavelengths, transparent_n);
        models["cauchy"] = {{"A", cauchy.A}, {"B", cauchy.B}, {"C", cauchy.C},
                            {"rms_deviation", cauchy.rms}, {"largest_deviation", cauchy.largest}};
        try {
          const auto sellmeier = ellidft::ellipsometry::fit_sellmeier(transparent_wavelengths, transparent_n);
          models["sellmeier"] = {{"B", {sellmeier.B}}, {"C_um2", {sellmeier.C_um2}}, {"converged", sellmeier.converged},
                                 {"rms_deviation", sellmeier.rms}, {"largest_deviation", sellmeier.largest}};
        } catch (const std::exception& e) {
          models["sellmeier"] = {{"error", e.what()}};
        }
      }
      out["dispersion_models"] = models;
    }
    std::cout << out.dump(2) << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
