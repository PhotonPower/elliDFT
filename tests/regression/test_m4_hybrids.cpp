// Milestone M4: erf-attenuated integrals and hybrid functionals against the
// reference code. Acceptance criteria:
//   - erf-attenuated integrals within 1e-10, plus the limits omega -> infinity
//     (Coulomb integrals) and omega -> 0 (zero);
//   - PBE0, B3LYP and CAM-B3LYP total energies within 5e-6 Hartree.
//
// Reference: tests/reference/data/m4/*.json (PySCF, tests/reference/gen_m4.py).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "df/df_integrals.hpp"
#include "integrals/shell_pair.hpp"
#include "integrals/two_electron.hpp"
#include "reference_data.hpp"
#include "scf/rhf.hpp"

using namespace ellidft;
using Catch::Matchers::WithinAbs;
using nlohmann::json;

namespace {

// Unique four-centre integrals in the packing of the reference
// (ij = i(i+1)/2 + j, position ij(ij+1)/2 + kl), each quartet computed once.
std::vector<double> packed_eri(const BasisSet& basis, double omega) {
  const auto n = static_cast<std::size_t>(basis.n_functions());
  const std::size_t n_pair = n * (n + 1) / 2;
  std::vector<double> packed(n_pair * (n_pair + 1) / 2, 0.0);
  const auto& shells = basis.shells();
  EriEngine engine(basis.max_l());
  for (std::size_t a = 0; a < shells.size(); ++a)
    for (std::size_t b = 0; b <= a; ++b) {
      const ShellPair bra = make_shell_pair(shells[a], shells[b]);
      for (std::size_t c = 0; c < shells.size(); ++c)
        for (std::size_t d = 0; d <= c; ++d) {
          const ShellPair ket = make_shell_pair(shells[c], shells[d]);
          const std::vector<double>& block = engine.compute(bra, ket, omega);
          std::size_t index = 0;
          for (int fa = 0; fa < bra.nfa; ++fa)
            for (int fb = 0; fb < bra.nfb; ++fb)
              for (int fc = 0; fc < ket.nfa; ++fc)
                for (int fd = 0; fd < ket.nfb; ++fd, ++index) {
                  const auto i = static_cast<std::size_t>(basis.first_function(a) + fa);
                  const auto j = static_cast<std::size_t>(basis.first_function(b) + fb);
                  const auto k = static_cast<std::size_t>(basis.first_function(c) + fc);
                  const auto l = static_cast<std::size_t>(basis.first_function(d) + fd);
                  if (j > i || l > k) continue;
                  const std::size_t ij = i * (i + 1) / 2 + j, kl = k * (k + 1) / 2 + l;
                  if (kl > ij) continue;
                  packed[ij * (ij + 1) / 2 + kl] = block[index];
                }
        }
    }
  return packed;
}

}  // namespace

TEST_CASE("M4: erf-attenuated integrals equal the reference", "[m4][regression]") {
  const auto files = testing::reference_files("m4", "erf_ints_");
  REQUIRE(files.size() == 3);
  for (const std::string& name : files) {
    DYNAMIC_SECTION(name) {
      const json record = testing::load_reference("m4/" + name);
      const double omega = record.at("omega").get<double>();
      const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), name);
      const BasisSet basis = testing::basis_from_reference(record, molecule);
      const BasisSet aux =
          BasisSet::build(molecule, BasisLibrary::load(record.at("auxiliary_basis").get<std::string>()), true);
      REQUIRE(basis.n_functions() == record.at("n_functions").get<int>());

      // Four centres.
      const std::vector<double> four = packed_eri(basis, omega);
      const auto reference = record.at("eri_packed").get<std::vector<double>>();
      REQUIRE(four.size() == reference.size());
      double dev4 = 0.0;
      for (std::size_t k = 0; k < four.size(); ++k) dev4 = std::max(dev4, std::abs(four[k] - reference[k]));

      // Two centres.
      const Eigen::MatrixXd V = two_center_coulomb(aux, omega);
      const double dev2 = (V - testing::unpack_symmetric(record.at("two_center"), V.rows())).cwiseAbs().maxCoeff();

      // Three centres (stored for the small basis sets only).
      double dev3 = 0.0;
      if (!record.at("three_center").is_null()) {
        const Eigen::MatrixXd A = three_center_coulomb(basis, aux, omega);
        const json& ref3 = record.at("three_center");
        for (Eigen::Index r = 0; r < A.rows(); ++r)
          for (Eigen::Index c = 0; c < A.cols(); ++c)
            dev3 = std::max(dev3, std::abs(A(r, c) - ref3[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)].get<double>()));
      }
      INFO("four-centre " << dev4 << ", three-centre " << dev3 << ", two-centre " << dev2);
      CHECK(dev4 <= 1e-10);
      CHECK(dev3 <= 1e-10);
      CHECK(dev2 <= 1e-10);
      SUCCEED(name << ": four-centre " << dev4 << ", three-centre " << dev3 << ", two-centre " << dev2);
    }
  }
}

namespace {

struct Summary {
  int cases = 0;
  double energy = 0.0;
};

void compare_energies(bool slow_cases, Summary& summary) {
  const json reference = testing::load_reference("m4/hybrids.json");
  const json& settings = reference.at("meta").at("settings");
  Method method;
  method.grid = {settings.at("grid").at("n_radial").get<int>(), settings.at("grid").at("n_angular").get<int>()};
  ScfOptions options;
  options.energy_tolerance = 1e-10;
  options.error_tolerance = 1e-7;

  // The exchange coefficients of both programs must be the same.
  for (const auto& [label, c] : settings.at("exchange_coefficients_in_pyscf").items()) {
    const Functional f = Functional::by_name(label);
    INFO(label);
    CHECK_THAT(f.exact_exchange, WithinAbs(c.at("exact_exchange").get<double>(), 1e-12));
    CHECK_THAT(f.long_range_exchange, WithinAbs(c.at("long_range_exchange").get<double>(), 1e-12));
    CHECK_THAT(f.omega, WithinAbs(c.at("omega").get<double>(), 1e-12));
  }

  for (const json& record : reference.at("results")) {
    const std::string label = record.at("molecule").get<std::string>() + " / " +
                              record.at("functional").get<std::string>();
    const Molecule molecule = testing::molecule_from_reference(record, label);
    const BasisSet basis = BasisSet::build(molecule, BasisLibrary::load(record.at("basis").get<std::string>()),
                                           record.at("spherical").get<bool>());
    const bool quick = molecule.n_atoms() <= 2 && basis.n_primitives() <= 12;
    if (quick == slow_cases) continue;
    INFO(label);
    method.functional = Functional::by_name(record.at("functional").get<std::string>());
    const ScfResult scf = run_scf(molecule, basis, method, options);
    CHECK(scf.converged);
    const double deviation = std::abs(scf.energy_total - record.at("energy_total").get<double>());
    CHECK(deviation <= 5e-6);  // acceptance criterion
    const auto dipole = record.at("dipole_au").get<std::vector<double>>();
    for (std::size_t k = 0; k < 3; ++k) CHECK(std::abs(scf.dipole[k] - dipole[k]) <= 1e-5);
    ++summary.cases;
    summary.energy = std::max(summary.energy, deviation);
    std::printf("HYB %-8s %-10s %2d iterations  E = %.10f  |dE| = %.2e\n",
                record.at("molecule").get<std::string>().c_str(), record.at("functional").get<std::string>().c_str(),
                scf.iterations, scf.energy_total, deviation);
    std::fflush(stdout);
  }
}

}  // namespace

TEST_CASE("M4: PBE0, B3LYP and CAM-B3LYP energies equal the reference (atoms and H2)", "[m4][regression]") {
  Summary summary;
  compare_energies(false, summary);
  CHECK(summary.cases == 9);
  SUCCEED("largest deviation over " << summary.cases << " cases: " << summary.energy << " Hartree");
}

TEST_CASE("M4: PBE0, B3LYP and CAM-B3LYP energies equal the reference (all other cases)",
          "[m4][regression][slow]") {
  Summary summary;
  compare_energies(true, summary);
  CHECK(summary.cases == 27);
  SUCCEED("largest deviation over " << summary.cases << " cases: " << summary.energy << " Hartree");
}
