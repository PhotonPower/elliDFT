// ellidft-ellipsometry: evaluation of ellipsometry measurements with a layer
// model (homogeneous, isotropic layers on a substrate).
//
// Usage: ellidft-ellipsometry simulate <model.json> [--table <file>]
//        ellidft-ellipsometry fit      <model.json> [--table <file>]
//
// The result is written as JSON to standard output, a short protocol to
// standard error; --table writes measured and calculated Psi and Delta as
// text. File names inside the model file are relative to its directory.
//
// Model file (docs/Ellipsometrie.md describes it in full):
//
//   {
//     "ambient":   {"model": "constant", "n": 1.0},
//     "layers": [                                   // from top to bottom
//       {"name": "film", "thickness_nm": {"value": 100, "fit": true, "scan": [50, 500, 10]},
//        "material": {"model": "cauchy", "A": {"value": 1.55, "fit": true}, "B": 0.007}},
//       {"name": "oxide", "thickness_nm": 2.0, "material": {"model": "cauchy", "A": 1.45, "B": 0.0036}}
//     ],
//     "substrate": {"model": "table", "file": "substrate.txt"},
//     "data":      {"file": "measurement.txt", "delta_offset_deg": 0},
//     "residuals": "ncs",
//     "simulate":  {"wavelengths_nm": {"from": 400, "to": 800, "step": 10}, "angles_deg": [65, 70]}
//   }
//
// Every number of a model may be given plainly (fixed) or as an object with
// "value", "fit", "min", "max" and "scan": [from, to, number of start values].

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "ellipsometry/analysis.hpp"

using namespace ellidft::ellipsometry;
using nlohmann::json;
using nlohmann::ordered_json;

namespace {

struct Model {
  Stack stack;
  std::vector<std::string> layer_names;
  std::vector<FitParameter> parameters;
  std::filesystem::path directory;
  // materials with a name ("materials" of the model file): definitions, and
  // the objects built so far. A material referred to in several places is one
  // object with one set of adjustable numbers.
  json definitions = json::object();
  std::map<std::string, std::shared_ptr<Material>> named;
  std::vector<std::string> in_progress;
};

// A number of the model: plain, or {"value", "fit", "min", "max", "scan"}.
// Returns the value; registers a fit parameter through `make` if requested.
double read_quantity(const json& node, const std::string& name, Model& model,
                     const std::function<FitParameter()>& make, double default_lower, double default_upper) {
  if (node.is_number()) return node.get<double>();
  if (!node.is_object() || !node.contains("value"))
    throw std::runtime_error(name + ": expected a number or an object with \"value\"");
  for (const auto& [key, unused] : node.items())
    if (key != "value" && key != "fit" && key != "min" && key != "max" && key != "scan")
      throw std::runtime_error(name + ": unknown key \"" + key + "\"");
  const double value = node.at("value").get<double>();
  if (node.value("fit", false)) {
    FitParameter p = make();
    p.name = name;
    p.lower = node.value("min", default_lower);
    p.upper = node.value("max", default_upper);
    if (node.contains("scan")) {
      const json& scan = node.at("scan");
      if (!scan.is_array() || scan.size() != 3) throw std::runtime_error(name + ": \"scan\" expects [from, to, number]");
      const double from = scan[0].get<double>(), to = scan[1].get<double>();
      const int number = scan[2].get<int>();
      if (number < 2) throw std::runtime_error(name + ": \"scan\" needs at least 2 start values");
      for (int i = 0; i < number; ++i) p.alternative_starts.push_back(from + (to - from) * i / (number - 1));
    }
    model.parameters.push_back(std::move(p));
  } else if (node.contains("scan") || node.contains("min") || node.contains("max")) {
    throw std::runtime_error(name + ": \"min\", \"max\" and \"scan\" need \"fit\": true");
  }
  return value;
}

// Table with wavelength_nm and two more columns (n k, or Re Im of alpha').
// `column_a` and `column_b` count from 1; column 1 is the wavelength.
void read_table(const std::filesystem::path& path, std::size_t column_a, std::size_t column_b, std::vector<double>& wavelengths,
                std::vector<Complex>& values) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path.string());
  std::string line;
  int line_number = 0;
  while (std::getline(in, line)) {
    ++line_number;
    const auto first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') continue;
    for (char& c : line)
      if (c == ',' || c == ';') c = ' ';
    std::istringstream fields(line);
    std::vector<double> numbers;
    for (double x = 0.0; fields >> x;) numbers.push_back(x);
    if (numbers.size() < column_a || numbers.size() < 2)
      throw std::runtime_error(path.string() + ", line " + std::to_string(line_number) + ": expected wavelength_nm n [k], or the columns asked for");
    wavelengths.push_back(numbers[0]);
    values.emplace_back(numbers[column_a - 1], numbers.size() >= column_b ? numbers[column_b - 1] : 0.0);
  }
}

std::shared_ptr<Material> read_material(const json& node, const std::string& name, Model& model) {
  if (node.is_string()) {   // reference to a named material
    const std::string reference = node.get<std::string>();
    if (const auto found = model.named.find(reference); found != model.named.end()) return found->second;
    if (!model.definitions.contains(reference)) throw std::runtime_error(name + ": no material named \"" + reference + "\"");
    if (std::find(model.in_progress.begin(), model.in_progress.end(), reference) != model.in_progress.end())
      throw std::runtime_error("materials: \"" + reference + "\" refers to itself");
    model.in_progress.push_back(reference);
    const std::shared_ptr<Material> material = read_material(model.definitions.at(reference), reference, model);
    model.in_progress.pop_back();
    model.named[reference] = material;
    return material;
  }
  if (!node.is_object() || !node.contains("model")) throw std::runtime_error(name + ": a material needs \"model\"");
  const std::string kind = node.at("model").get<std::string>();
  auto check_keys = [&](std::initializer_list<const char*> allowed) {
    for (const auto& [key, unused] : node.items()) {
      bool known = key == "model";
      for (const char* a : allowed) known = known || key == a;
      if (!known) throw std::runtime_error(name + " (" + kind + "): unknown key \"" + key + "\"");
    }
  };
  const double inf = std::numeric_limits<double>::infinity();
  // Reads the adjustable numbers of a material that already exists with
  // default values, then sets them.
  auto quantities = [&](const std::shared_ptr<Material>& material, const std::vector<std::pair<std::string, double>>& defaults,
                        const std::vector<std::pair<double, double>>& bounds) {
    for (std::size_t i = 0; i < defaults.size(); ++i) {
      const std::string& key = defaults[i].first;
      double value = defaults[i].second;
      if (node.contains(key))
        value = read_quantity(node.at(key), name + "." + key, model,
                              [material, i] { return material_parameter(material, i, ""); }, bounds[i].first, bounds[i].second);
      else if (std::isnan(value))
        throw std::runtime_error(name + " (" + kind + "): \"" + key + "\" is missing");
      material->set_parameter(i, value);
    }
  };
  const double required = std::numeric_limits<double>::quiet_NaN();

  if (kind == "constant") {
    check_keys({"n", "k"});
    auto material = std::make_shared<ConstantMaterial>(1.0, 0.0);
    quantities(material, {{"n", required}, {"k", 0.0}}, {{0.0, inf}, {0.0, inf}});
    return material;
  }
  if (kind == "cauchy") {
    check_keys({"A", "B", "C"});
    auto material = std::make_shared<CauchyMaterial>(1.0, 0.0, 0.0);
    quantities(material, {{"A", required}, {"B", 0.0}, {"C", 0.0}}, {{1.0, inf}, {-inf, inf}, {-inf, inf}});
    return material;
  }
  if (kind == "sellmeier") {
    check_keys({"B", "C_um2"});
    const json& B = node.at("B");
    const json& C = node.at("C_um2");
    if (!B.is_array() || !C.is_array() || B.size() != C.size() || B.empty())
      throw std::runtime_error(name + " (sellmeier): \"B\" and \"C_um2\" must be lists of equal length");
    auto material = std::make_shared<SellmeierMaterial>(std::vector<double>(B.size(), 0.0), std::vector<double>(B.size(), 0.0));
    for (std::size_t i = 0; i < B.size(); ++i) {
      const std::string term = std::to_string(i + 1);
      material->set_parameter(2 * i, read_quantity(B[i], name + ".B" + term, model,
                                                   [material, i] { return material_parameter(material, 2 * i, ""); }, 0.0, inf));
      material->set_parameter(2 * i + 1, read_quantity(C[i], name + ".C" + term, model,
                                                       [material, i] { return material_parameter(material, 2 * i + 1, ""); }, 0.0, inf));
    }
    return material;
  }
  if (kind == "table") {
    check_keys({"file", "columns"});
    std::size_t column_a = 2, column_b = 3;
    if (node.contains("columns")) {
      const std::vector<int> columns = node.at("columns").get<std::vector<int>>();
      if (columns.size() != 2 || columns[0] < 2 || columns[1] < 2)
        throw std::runtime_error(name + " (table): \"columns\" expects two column numbers from 2 on (n and k)");
      column_a = static_cast<std::size_t>(columns[0]);
      column_b = static_cast<std::size_t>(columns[1]);
    }
    std::vector<double> wavelengths;
    std::vector<Complex> indices;
    read_table(model.directory / node.at("file").get<std::string>(), column_a, column_b, wavelengths, indices);
    return std::make_shared<TabulatedMaterial>(std::move(wavelengths), std::move(indices));
  }
  if (kind == "ellidft") {
    // Output of ellidft-optics: polarisability volume per repeat unit at the
    // computed wavelengths; the density is taken from that file unless given.
    check_keys({"file", "density_g_cm3"});
    const std::filesystem::path path = model.directory / node.at("file").get<std::string>();
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read " + path.string());
    const json optics = json::parse(in);
    std::vector<std::pair<double, Complex>> rows;
    for (const json& row : optics.at("optical_constants")) {
      const json& alpha = row.at("increment_mean_A3");
      rows.emplace_back(row.at("wavelength_nm").get<double>(), Complex(alpha.at(0).get<double>(), alpha.at(1).get<double>()));
    }
    std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<double> wavelengths;
    std::vector<Complex> alpha;
    for (const auto& [w, a] : rows) {
      wavelengths.push_back(w);
      alpha.push_back(a);
    }
    const json& info = optics.at("material");
    auto material = std::make_shared<LorentzLorenzMaterial>(std::move(wavelengths), std::move(alpha),
                                                            info.at("molar_mass_repeat_unit_g_mol").get<double>(),
                                                            info.at("density_g_cm3").get<double>());
    quantities(material, {{"density_g_cm3", info.at("density_g_cm3").get<double>()}}, {{0.0, inf}});
    return material;
  }
  if (kind == "lorentz") {
    check_keys({"eps_inf", "oscillators"});
    const json& list = node.at("oscillators");
    if (!list.is_array() || list.empty()) throw std::runtime_error(name + " (lorentz): \"oscillators\" must be a list");
    auto material = std::make_shared<LorentzMaterial>(1.0, std::vector<LorentzMaterial::Oscillator>(list.size(), {0.0, 1.0, 0.1}));
    material->set_parameter(0, node.contains("eps_inf")
                                   ? read_quantity(node.at("eps_inf"), name + ".eps_inf", model,
                                                   [material] { return material_parameter(material, 0, ""); }, 1.0, inf)
                                   : 1.0);
    for (std::size_t j = 0; j < list.size(); ++j) {
      for (const auto& [key, unused] : list[j].items())
        if (key != "strength_ev2" && key != "energy_ev" && key != "width_ev")
          throw std::runtime_error(name + " (lorentz): unknown key \"" + key + "\" in an oscillator");
      const char* keys[3] = {"strength_ev2", "energy_ev", "width_ev"};
      for (std::size_t k = 0; k < 3; ++k) {
        const std::size_t index = 1 + 3 * j + k;
        material->set_parameter(index, read_quantity(list[j].at(keys[k]), name + "." + keys[k] + std::to_string(j + 1), model,
                                                     [material, index] { return material_parameter(material, index, ""); }, 0.0, inf));
      }
    }
    return material;
  }
  if (kind == "tauc_lorentz") {
    check_keys({"eps_inf", "A_ev", "E0_ev", "C_ev", "Eg_ev"});
    auto material = std::make_shared<TaucLorentzMaterial>(1.0, 0.0, 1.0, 1.0, 0.0);
    quantities(material, {{"eps_inf", 1.0}, {"A_ev", required}, {"E0_ev", required}, {"C_ev", required}, {"Eg_ev", required}},
               {{0.0, inf}, {0.0, inf}, {1e-3, inf}, {1e-3, inf}, {0.0, inf}});
    return material;
  }
  if (kind == "bruggeman") {
    check_keys({"a", "b", "fraction_b"});
    const std::shared_ptr<Material> a = read_material(node.at("a"), name + ".a", model);
    const std::shared_ptr<Material> b = read_material(node.at("b"), name + ".b", model);
    auto material = std::make_shared<BruggemanMaterial>(a, b, 0.5);
    quantities(material, {{"fraction_b", required}}, {{0.0, 1.0}});
    return material;
  }
  throw std::runtime_error(name + ": unknown material model \"" + kind +
                           "\" (constant, cauchy, sellmeier, lorentz, tauc_lorentz, bruggeman, table, ellidft)");
}

std::unique_ptr<Model> read_model(const json& input, const std::filesystem::path& directory) {
  for (const auto& [key, unused] : input.items())
    if (key != "materials" && key != "substrate_thickness_mm" && key != "backside" &&
        key != "ambient" && key != "layers" && key != "substrate" && key != "data" && key != "residuals" &&
        key != "simulate" && key != "comment")
      throw std::runtime_error("model file: unknown key \"" + key + "\"");
  auto model = std::make_unique<Model>();   // fit parameters keep references into the stack
  model->directory = directory;
  if (input.contains("materials")) {
    if (!input.at("materials").is_object()) throw std::runtime_error("model file: \"materials\" must be an object of named materials");
    model->definitions = input.at("materials");
  }
  model->stack.ambient = input.contains("ambient") ? read_material(input.at("ambient"), "ambient", *model)
                                                    : std::make_shared<ConstantMaterial>(1.0);
  if (!input.contains("substrate")) throw std::runtime_error("model file: \"substrate\" is missing");
  const json layers = input.value("layers", json::array());
  model->stack.layers.resize(layers.size());
  for (std::size_t j = 0; j < layers.size(); ++j) {
    const json& layer = layers[j];
    const std::string name = layer.value("name", "layer" + std::to_string(j + 1));
    for (const auto& [key, unused] : layer.items())
      if (key != "name" && key != "thickness_nm" && key != "material" && key != "material_extraordinary")
        throw std::runtime_error(name + ": unknown key \"" + key + "\"");
    model->layer_names.push_back(name);
    model->stack.layers[j].material = read_material(layer.at("material"), name, *model);
    if (layer.contains("material_extraordinary"))
      model->stack.layers[j].extraordinary = read_material(layer.at("material_extraordinary"), name + ".extraordinary", *model);
    Stack& stack = model->stack;
    model->stack.layers[j].thickness_nm =
        read_quantity(layer.at("thickness_nm"), name + ".thickness_nm", *model,
                      [&stack, j] { return thickness_parameter(stack, j, ""); }, 0.0, std::numeric_limits<double>::infinity());
    if (model->stack.layers[j].thickness_nm < 0.0) throw std::runtime_error(name + ": negative thickness");
  }
  model->stack.substrate = read_material(input.at("substrate"), "substrate", *model);
  if (input.contains("substrate_thickness_mm")) {
    const double mm = input.at("substrate_thickness_mm").get<double>();
    if (!(mm > 0.0)) throw std::runtime_error("model file: \"substrate_thickness_mm\" must be positive");
    model->stack.substrate_thickness_nm = mm * 1e6;
  }
  if (input.contains("backside")) {
    if (!input.contains("substrate_thickness_mm"))
      throw std::runtime_error("model file: \"backside\" needs \"substrate_thickness_mm\"");
    model->stack.backside = read_material(input.at("backside"), "backside", *model);
  }
  return model;
}

std::vector<double> read_range(const json& node, const std::string& name) {
  std::vector<double> values;
  if (node.is_array()) {
    values = node.get<std::vector<double>>();
  } else if (node.is_number()) {
    values.push_back(node.get<double>());
  } else {
    const double from = node.at("from").get<double>(), to = node.at("to").get<double>(), step = node.at("step").get<double>();
    if (!(step > 0.0) || to < from) throw std::runtime_error(name + ": need from <= to and step > 0");
    for (int i = 0; from + i * step <= to + 1e-9 * step; ++i) values.push_back(from + i * step);
  }
  if (values.empty()) throw std::runtime_error(name + ": no values");
  return values;
}

ordered_json layers_json(const Model& model) {
  ordered_json layers = ordered_json::array();
  for (std::size_t j = 0; j < model.stack.layers.size(); ++j) {
    ordered_json layer = {{"name", model.layer_names[j]}, {"thickness_nm", model.stack.layers[j].thickness_nm}};
    const Material& material = *model.stack.layers[j].material;
    const auto names = material.parameter_names();
    for (std::size_t i = 0; i < names.size(); ++i) layer["material"][names[i]] = material.parameter(i);
    layers.push_back(layer);
  }
  return layers;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::string command, model_file, table_file;
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--table") {
        if (i + 1 >= argc) throw std::invalid_argument("option --table needs a value");
        table_file = argv[++i];
      } else if (!arg.empty() && arg[0] == '-') {
        throw std::invalid_argument("unknown option " + arg);
      } else if (command.empty()) {
        command = arg;
      } else if (model_file.empty()) {
        model_file = arg;
      } else {
        throw std::invalid_argument("unexpected argument " + arg);
      }
    }
    if ((command != "simulate" && command != "fit") || model_file.empty()) {
      std::cerr << "usage: ellidft-ellipsometry simulate|fit <model.json> [--table <file>]\n";
      return EXIT_FAILURE;
    }

    std::ifstream in(model_file);
    if (!in) throw std::runtime_error("cannot read " + model_file);
    const json input = json::parse(in, nullptr, true, true);   // comments allowed
    const auto model = read_model(input, std::filesystem::absolute(model_file).parent_path());

    ordered_json out;
    out["program"] = "ellidft-ellipsometry";
    out["command"] = command;
    out["conventions"] = {{"rho", "r_p / r_s = tan(Psi) exp(i Delta)"},
                          {"delta", "180 degrees for a bare transparent substrate below the Brewster angle"},
                          {"refractive_index", "N = n + i k, k >= 0 for absorption"},
                          {"model", "homogeneous layers, isotropic or uniaxial with the axis along the normal; sharp "
                                    "interfaces; back face of a thick substrate added incoherently if a thickness is given"}};

    std::vector<Measurement> points;
    if (command == "simulate") {
      if (!model->parameters.empty()) std::cerr << "note: \"fit\" entries are ignored by simulate\n";
      if (!input.contains("simulate")) throw std::runtime_error("model file: \"simulate\" is missing");
      const std::vector<double> wavelengths = read_range(input.at("simulate").at("wavelengths_nm"), "simulate.wavelengths_nm");
      const std::vector<double> angles = read_range(input.at("simulate").at("angles_deg"), "simulate.angles_deg");
      for (double angle : angles)
        for (double wavelength : wavelengths) {
          Measurement m;
          m.wavelength_nm = wavelength;
          m.angle_deg = angle;
          points.push_back(m);
        }
      out["layers"] = layers_json(*model);
    } else {
      if (!input.contains("data")) throw std::runtime_error("model file: \"data\" is missing");
      const json& data = input.at("data");
      const std::filesystem::path path = model->directory / data.at("file").get<std::string>();
      points = read_measurements_file(path.string(), data.value("delta_offset_deg", 0.0));
      const std::string residual_name = input.value("residuals", "ncs");
      if (residual_name != "ncs" && residual_name != "psi_delta")
        throw std::runtime_error("model file: \"residuals\" expects \"ncs\" or \"psi_delta\"");
      const Residuals kind = residual_name == "ncs" ? Residuals::NCS : Residuals::PsiDelta;
      if (model->parameters.empty()) throw std::runtime_error("model file: no quantity has \"fit\": true");

      const Deviations before = deviations(model->stack, points);
      const StackFit result = fit_stack(model->stack, model->parameters, points, kind);

      std::cerr << "fit of " << model->parameters.size() << " parameters to " << points.size() << " points, "
                << result.n_starts << " start(s), residuals " << residual_name << "\n"
                << "  " << (result.fit.converged ? "converged" : "NOT CONVERGED") << " after " << result.fit.iterations
                << " iterations (" << result.fit.message << ")\n";
      ordered_json parameters = ordered_json::array();
      for (std::size_t j = 0; j < result.names.size(); ++j) {
        const auto J = static_cast<Eigen::Index>(j);
        const double value = result.fit.parameters[J];
        const FitParameter& p = model->parameters[j];
        const bool at_bound = value <= p.lower || value >= p.upper;
        ordered_json entry = {{"name", result.names[j]}, {"value", value}, {"at_bound", at_bound}};
        if (result.fit.covariance_valid) entry["standard_error"] = result.fit.standard_errors[J];
        else entry["standard_error"] = nullptr;
        parameters.push_back(entry);
        std::cerr << "  " << std::left << std::setw(28) << result.names[j] << std::right << std::setprecision(8)
                  << std::setw(16) << value;
        if (result.fit.covariance_valid) std::cerr << "  +- " << std::setprecision(3) << result.fit.standard_errors[J];
        if (at_bound) std::cerr << "  (at its bound)";
        std::cerr << "\n";
      }
      std::cerr << std::setprecision(4) << "  rms deviation: Psi " << result.deviation.rms_psi_deg << " deg, Delta "
                << result.deviation.rms_delta_deg << " deg (before the fit: " << before.rms_psi_deg << ", "
                << before.rms_delta_deg << ")\n";
      if (!result.fit.covariance_valid)
        std::cerr << "  WARNING: the parameters are not independent (or one has no influence); no error estimates\n";

      ordered_json fit = {{"converged", result.fit.converged},
                          {"message", result.fit.message},
                          {"iterations", result.fit.iterations},
                          {"model_evaluations", result.fit.evaluations},
                          {"starts", result.n_starts},
                          {"residuals", residual_name},
                          {"n_points", points.size()},
                          {"sum_of_squares", result.fit.sum_of_squares},
                          {"rms_residual_per_degree_of_freedom", result.rms_residual},
                          {"rms_psi_deg", result.deviation.rms_psi_deg},
                          {"rms_delta_deg", result.deviation.rms_delta_deg},
                          {"max_psi_deg", result.deviation.max_psi_deg},
                          {"max_delta_deg", result.deviation.max_delta_deg},
                          {"rms_ncs", result.deviation.rms_ncs},
                          {"parameters", parameters},
                          {"errors", "standard errors from the curvature at the minimum, scatter estimated from the "
                                     "fit itself; they do not include errors of the model or of fixed quantities"}};
      if (result.fit.covariance_valid) {
        ordered_json correlation = ordered_json::array();
        for (Eigen::Index a = 0; a < result.fit.correlation.rows(); ++a) {
          std::vector<double> row(static_cast<std::size_t>(result.fit.correlation.cols()));
          for (Eigen::Index b = 0; b < result.fit.correlation.cols(); ++b) row[static_cast<std::size_t>(b)] = result.fit.correlation(a, b);
          correlation.push_back(row);
        }
        fit["correlation"] = correlation;
      }
      out["fit"] = fit;
      out["layers"] = layers_json(*model);
    }

    const std::vector<PsiDelta> calculated = simulate(model->stack, points);
    ordered_json rows = ordered_json::array();
    std::ofstream table;
    if (!table_file.empty()) {
      table.open(table_file);
      if (!table) throw std::runtime_error("cannot write " + table_file);
      table << (command == "fit" ? "# wavelength_nm angle_deg psi_meas delta_meas psi_calc delta_calc\n"
                                 : "# wavelength_nm angle_deg psi_deg delta_deg\n");
      table << std::fixed;
    }
    for (std::size_t i = 0; i < points.size(); ++i) {
      ordered_json row = {{"wavelength_nm", points[i].wavelength_nm}, {"angle_deg", points[i].angle_deg}};
      if (command == "fit") {
        row["psi_measured_deg"] = points[i].psi_deg;
        row["delta_measured_deg"] = points[i].delta_deg;
      }
      row["psi_deg"] = calculated[i].psi_deg;
      row["delta_deg"] = calculated[i].delta_deg;
      rows.push_back(row);
      if (table.is_open()) {
        table << std::setprecision(3) << points[i].wavelength_nm << " " << points[i].angle_deg << std::setprecision(5);
        if (command == "fit") table << " " << points[i].psi_deg << " " << points[i].delta_deg;
        table << " " << calculated[i].psi_deg << " " << calculated[i].delta_deg << "\n";
      }
    }
    out["points"] = rows;

    // Optical constants of the layers at the wavelengths used.
    std::vector<double> wavelengths;
    for (const Measurement& m : points)
      if (std::find(wavelengths.begin(), wavelengths.end(), m.wavelength_nm) == wavelengths.end()) wavelengths.push_back(m.wavelength_nm);
    ordered_json constants = ordered_json::array();
    for (std::size_t j = 0; j < model->stack.layers.size(); ++j) {
      ordered_json entry = {{"name", model->layer_names[j]}, {"n_k", ordered_json::array()}};
      for (double wavelength : wavelengths) {
        const Complex N = model->stack.layers[j].material->index(wavelength);
        entry["n_k"].push_back({wavelength, N.real(), N.imag()});
      }
      constants.push_back(entry);
    }
    out["optical_constants"] = constants;

    std::cout << out.dump(2) << "\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return EXIT_FAILURE;
  }
}
