// Milestone M5: linear response against the reference code.
// Acceptance criteria, for PBE, PBE0 and CAM-B3LYP:
//   - static polarisability within a relative 1e-4;
//   - the ten lowest singlet excitations, full response and Tamm-Dancoff:
//     energies within 1e-5 Hartree, oscillator strengths within 1e-4.
//
// Reference: tests/reference/data/m5/response.json (PySCF, tests/reference/gen_m5.py).
//
// Every test run covers hydrogen fluoride with PBE0 (a few seconds). The other
// cases are tagged [slow] and meant for a server.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "reference_data.hpp"
#include "response/excitations.hpp"
#include "response/polarizability.hpp"
#include "response/response_operator.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using nlohmann::json;

namespace {

struct Summary {
  int cases = 0;
  double energy = 0.0, strength = 0.0, polarizability = 0.0;
};

// Oscillator strengths of degenerate states can be distributed differently
// among the partners; they are compared as sums over groups of states whose
// reference energies agree within 1e-5 Hartree.
void compare_states(const ExcitationResult& result, const json& reference, Summary& summary) {
  const auto energies = reference.at("energies").get<std::vector<double>>();
  const auto strengths = reference.at("oscillator_strengths").get<std::vector<double>>();
  REQUIRE(result.states.size() == energies.size());
  std::size_t start = 0;
  while (start < energies.size()) {
    std::size_t end = start + 1;
    while (end < energies.size() && energies[end] - energies[end - 1] < 1e-5) ++end;
    const bool group_complete = end < energies.size();  // the last group may be cut by the number of states
    double f = 0.0, f_reference = 0.0;
    for (std::size_t s = start; s < end; ++s) {
      INFO("state " << s);
      const double deviation = std::abs(result.states[s].energy - energies[s]);
      CHECK(deviation <= 1e-5);  // acceptance criterion
      summary.energy = std::max(summary.energy, deviation);
      f += result.states[s].oscillator_strength;
      f_reference += strengths[s];
    }
    if (group_complete || end - start == 1) {
      INFO("states " << start << " to " << end - 1);
      CHECK(std::abs(f - f_reference) <= 1e-4);  // acceptance criterion
      summary.strength = std::max(summary.strength, std::abs(f - f_reference));
    }
    start = end;
  }
}

void compare_response(bool slow_cases, Summary& summary) {
  const json reference = testing::load_reference("m5/response.json");
  const json& settings = reference.at("meta").at("settings");
  Method method;
  method.grid = {settings.at("grid").at("n_radial").get<int>(), settings.at("grid").at("n_angular").get<int>()};
  ScfOptions scf_options;
  scf_options.energy_tolerance = 1e-11;
  scf_options.error_tolerance = 1e-9;

  for (const json& record : reference.at("results")) {
    const std::string molecule_name = record.at("molecule").get<std::string>();
    const std::string functional_name = record.at("functional").get<std::string>();
    const bool quick = molecule_name == "hf" && functional_name == "PBE0";
    if (quick == slow_cases) continue;
    const std::string label = molecule_name + " / " + functional_name;
    INFO(label);
    const Molecule molecule = testing::molecule_from_reference(record, label);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(record.at("basis").get<std::string>()),
                                           record.at("spherical").get<bool>());
    method.functional = Functional::by_name(functional_name);
    DirectJK jk(basis);
    const ScfResult scf = run_scf(molecule, basis, method, jk, scf_options);
    REQUIRE(scf.converged);
    CHECK(std::abs(scf.energy_total - record.at("energy_total").get<double>()) <= 5e-6);
    const ResponseOperator op(molecule, basis, method, scf, jk);

    ExcitationOptions options;
    options.n_states = settings.at("n_states").get<int>();
    options.tolerance = 1e-6;
    {
      INFO("full response");
      const ExcitationResult rpa = excitations(op, options);
      CHECK(rpa.converged);
      compare_states(rpa, record.at("rpa"), summary);
    }
    {
      INFO("Tamm-Dancoff");
      options.tamm_dancoff = true;
      const ExcitationResult tda = excitations(op, options);
      CHECK(tda.converged);
      compare_states(tda, record.at("tda"), summary);
    }

    if (!record.at("polarizability_static").is_null()) {
      const Polarizability alpha = polarizability(op);
      CHECK(alpha.converged);
      const auto ref = record.at("polarizability_static").get<std::vector<std::vector<double>>>();
      double largest = 0.0;
      for (const auto& row : ref)
        for (double v : row) largest = std::max(largest, std::abs(v));
      for (std::size_t k = 0; k < 3; ++k)
        for (std::size_t l = 0; l < 3; ++l) {
          INFO("alpha component " << k << l);
          // Relative to the component; components that vanish by symmetry are
          // compared on the scale of the largest component.
          const double scale = std::abs(ref[k][l]) > 1e-6 * largest ? std::abs(ref[k][l]) : largest;
          const double deviation = std::abs(alpha.real(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(l)) -
                                            ref[k][l]) / scale;
          CHECK(deviation <= 1e-4);  // acceptance criterion
          summary.polarizability = std::max(summary.polarizability, deviation);
        }
    }
    ++summary.cases;
    std::printf("RESPONSE %-8s %-10s done\n", molecule_name.c_str(), functional_name.c_str());
    std::fflush(stdout);
  }
}

}  // namespace

TEST_CASE("M5: excitations and polarisability equal the reference (HF molecule, PBE0)", "[m5][regression]") {
  Summary summary;
  compare_response(false, summary);
  CHECK(summary.cases == 1);
  SUCCEED("largest deviations: excitation energy " << summary.energy << " Hartree, oscillator strength "
                                                   << summary.strength << ", polarisability (relative) "
                                                   << summary.polarizability);
}

TEST_CASE("M5: excitations and polarisability equal the reference (all other cases)", "[m5][regression][slow]") {
  Summary summary;
  compare_response(true, summary);
  CHECK(summary.cases == 8);
  SUCCEED("largest deviations over " << summary.cases << " cases: excitation energy " << summary.energy
                                     << " Hartree, oscillator strength " << summary.strength
                                     << ", polarisability (relative) " << summary.polarizability);
}
