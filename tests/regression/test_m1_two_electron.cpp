// Milestone M1: two-electron integrals against the reference code.
// Acceptance criterion: full tensor for small molecules, largest deviation at
// most 1e-10, eightfold permutational symmetry fulfilled.
//
// Reference: tests/reference/data/m1/eri_*.json and highl_*.json
// (PySCF, tests/reference/gen_m1.py).

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "full_eri.hpp"
#include "integrals/shell_pair.hpp"
#include "integrals/two_electron.hpp"
#include "reference_data.hpp"

using namespace ellidft;
using nlohmann::json;

namespace {

constexpr double acceptance_threshold = 1e-10;

// Largest violation of (ij|kl) = (ji|kl) = (ij|lk) = (kl|ij) over all indices.
// These three generate the eightfold symmetry.
double symmetry_violation(const testing::FullEri& eri) {
  double worst = 0.0;
  for (int i = 0; i < eri.n; ++i)
    for (int j = 0; j < eri.n; ++j)
      for (int k = 0; k < eri.n; ++k)
        for (int l = 0; l < eri.n; ++l) {
          const double v = eri(i, j, k, l);
          worst = std::max({worst, std::abs(v - eri(j, i, k, l)), std::abs(v - eri(i, j, l, k)),
                            std::abs(v - eri(k, l, i, j))});
        }
  return worst;
}

// Largest deviation from the reference, which stores the unique integrals:
// ij = i(i+1)/2 + j for i >= j, position ij(ij+1)/2 + kl for ij >= kl.
double deviation_from_packed(const testing::FullEri& eri, const json& packed) {
  const auto n = static_cast<std::size_t>(eri.n);
  const std::size_t n_pair = n * (n + 1) / 2;
  REQUIRE(packed.size() == n_pair * (n_pair + 1) / 2);
  double worst = 0.0;
  for (int i = 0; i < eri.n; ++i)
    for (int j = 0; j <= i; ++j)
      for (int k = 0; k <= i; ++k)
        for (int l = 0; l <= k; ++l) {
          const std::size_t ij = static_cast<std::size_t>(i) * static_cast<std::size_t>(i + 1) / 2 +
                                 static_cast<std::size_t>(j);
          const std::size_t kl = static_cast<std::size_t>(k) * static_cast<std::size_t>(k + 1) / 2 +
                                 static_cast<std::size_t>(l);
          if (kl > ij) continue;
          worst = std::max(worst, std::abs(eri(i, j, k, l) - packed[ij * (ij + 1) / 2 + kl].get<double>()));
        }
  return worst;
}

}  // namespace

TEST_CASE("M1: full two-electron tensors equal the reference and are symmetric", "[m1][regression]") {
  const auto files = testing::reference_files("m1", "eri_");
  REQUIRE(files.size() == 11);
  double worst_deviation = 0.0, worst_symmetry = 0.0;
  for (const std::string& name : files) {
    const json record = testing::load_reference("m1/" + name);
    const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), "reference");
    const BasisSet basis = testing::basis_from_reference(record, molecule);
    REQUIRE(basis.n_functions() == record.at("n_functions").get<int>());

    const testing::FullEri eri = testing::full_eri_tensor(basis);
    const double deviation = deviation_from_packed(eri, record.at("eri_packed"));
    const double symmetry = symmetry_violation(eri);
    INFO(name << ": deviation " << deviation << ", symmetry violation " << symmetry);
    CHECK(deviation <= acceptance_threshold);
    CHECK(symmetry <= acceptance_threshold);
    worst_deviation = std::max(worst_deviation, deviation);
    worst_symmetry = std::max(worst_symmetry, symmetry);
  }
  SUCCEED("largest values over 11 cases: deviation " << worst_deviation << ", symmetry violation "
                                                    << worst_symmetry);
}

// Shell quartets with angular momenta up to l = 6 (total up to 19), on three
// centres without symmetry, spherical and cartesian.
TEST_CASE("M1: two-electron shell quartets up to l = 6 equal the reference", "[m1][regression]") {
  for (const std::string name : {"highl_sph.json", "highl_cart.json"}) {
    DYNAMIC_SECTION(name) {
      const json record = testing::load_reference("m1/" + name);
      const Molecule molecule = testing::molecule_from_reference(record.at("molecule"), "reference");
      const BasisSet basis = testing::basis_from_reference(record, molecule);
      REQUIRE(basis.n_shells() == record.at("n_shells").get<std::size_t>());
      const auto& shells = basis.shells();
      EriEngine engine(basis.max_l());

      double worst = 0.0;
      for (const json& quartet : record.at("eri_shell_quartets")) {
        const auto s = quartet.at("shells").get<std::vector<std::size_t>>();
        const auto shape = quartet.at("shape").get<std::vector<int>>();
        const ShellPair bra = make_shell_pair(shells[s[0]], shells[s[1]]);
        const ShellPair ket = make_shell_pair(shells[s[2]], shells[s[3]]);
        REQUIRE(bra.nfa == shape[0]);
        REQUIRE(bra.nfb == shape[1]);
        REQUIRE(ket.nfa == shape[2]);
        REQUIRE(ket.nfb == shape[3]);
        const std::vector<double>& block = engine.compute(bra, ket);
        const auto reference = quartet.at("values").get<std::vector<double>>();
        REQUIRE(block.size() == reference.size());
        double deviation = 0.0;
        for (std::size_t k = 0; k < block.size(); ++k)
          deviation = std::max(deviation, std::abs(block[k] - reference[k]));
        INFO("shells " << s[0] << " " << s[1] << " " << s[2] << " " << s[3] << " (l = " << bra.la << " "
                       << bra.lb << " " << ket.la << " " << ket.lb << "): deviation " << deviation);
        CHECK(deviation <= acceptance_threshold);
        worst = std::max(worst, deviation);
      }
      SUCCEED(name << ": largest deviation " << worst);
    }
  }
}
