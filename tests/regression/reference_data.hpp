#pragma once

// Access to the stored reference data (tests/reference/data), produced by the
// PySCF scripts in tests/reference. The tests run without PySCF.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <nlohmann/json.hpp>

#include "basis/basis_set.hpp"
#include "core/elements.hpp"
#include "core/molecule.hpp"

namespace ellidft::testing {

inline std::filesystem::path source_dir() { return ELLIDFT_SOURCE_DIR; }

inline nlohmann::json load_reference(const std::string& relative_path) {
  const auto path = source_dir() / "tests" / "reference" / "data" / relative_path;
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open reference file " + path.string());
  return nlohmann::json::parse(in);
}

// Builds a molecule from the "symbols" / "geometry_bohr" / "charge" entries of
// a reference record. Comparisons always use the geometry in Bohr.
inline Molecule molecule_from_reference(const nlohmann::json& record, const std::string& name) {
  Molecule mol;
  mol.name = name;
  mol.charge = record.at("charge").get<int>();
  const auto& symbols = record.at("symbols");
  const auto& geometry = record.at("geometry_bohr");
  if (symbols.size() != geometry.size())
    throw std::runtime_error("reference record: symbols and geometry differ in length");
  for (std::size_t i = 0; i < symbols.size(); ++i) {
    Atom atom;
    atom.Z = atomic_number(symbols[i].get<std::string>());
    if (atom.Z == 0)
      throw std::runtime_error("reference record: unknown element " + symbols[i].get<std::string>());
    for (std::size_t k = 0; k < 3; ++k) atom.r[k] = geometry[i].at(k).get<double>();
    mol.atoms.push_back(atom);
  }
  return mol;
}

// Symmetric matrix from its lower triangle stored row by row
// ((0,0) (1,0) (1,1) (2,0) ...), the packing of the reference files.
inline Eigen::MatrixXd unpack_symmetric(const nlohmann::json& packed, Eigen::Index n) {
  if (packed.size() != static_cast<std::size_t>(n * (n + 1) / 2))
    throw std::runtime_error("reference record: packed matrix has the wrong length");
  Eigen::MatrixXd m(n, n);
  std::size_t k = 0;
  for (Eigen::Index i = 0; i < n; ++i)
    for (Eigen::Index j = 0; j <= i; ++j) m(i, j) = m(j, i) = packed[k++].get<double>();
  return m;
}

// Basis of a reference record: by name from data/basis, or from the explicit
// shells of an artificial basis ("shells_per_atom").
inline BasisSet basis_from_reference(const nlohmann::json& record, const Molecule& molecule) {
  const nlohmann::json& basis = record.at("basis");
  const bool spherical = basis.at("spherical").get<bool>();
  if (!basis.contains("shells_per_atom"))
    return BasisSet::build(molecule, BasisLibrary::load(basis.at("name").get<std::string>()), spherical);
  std::vector<Shell> shells;
  const nlohmann::json& per_atom = basis.at("shells_per_atom");
  if (per_atom.size() != molecule.n_atoms())
    throw std::runtime_error("reference record: shells_per_atom does not match the molecule");
  for (std::size_t a = 0; a < per_atom.size(); ++a)
    for (const nlohmann::json& s : per_atom[a]) {
      Shell shell;
      shell.atom = static_cast<int>(a);
      shell.center = molecule.atoms[a].r;
      shell.l = s.at("l").get<int>();
      shell.exponents = s.at("exponents").get<std::vector<double>>();
      shell.coefficients = s.at("coefficients").get<std::vector<double>>();
      shells.push_back(std::move(shell));
    }
  return BasisSet::from_shells(basis.at("name").get<std::string>(), spherical, std::move(shells));
}

// Names of the reference files tests/reference/data/<directory>/<prefix>*.json, sorted.
inline std::vector<std::string> reference_files(const std::string& directory, const std::string& prefix) {
  std::vector<std::string> names;
  for (const auto& entry :
       std::filesystem::directory_iterator(source_dir() / "tests" / "reference" / "data" / directory)) {
    const std::string name = entry.path().filename().string();
    if (name.rfind(prefix, 0) == 0 && entry.path().extension() == ".json") names.push_back(name);
  }
  std::sort(names.begin(), names.end());
  return names;
}

}  // namespace ellidft::testing
