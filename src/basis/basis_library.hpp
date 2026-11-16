#pragma once

#include <filesystem>
#include <iosfwd>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace ellidft {

// Highest angular momentum the basis module accepts (s p d f g h i).
inline constexpr int max_angular_momentum = 6;

// Angular momentum for a shell letter (case-insensitive), -1 if unknown.
int angular_momentum_from_letter(char letter);

// Lower-case shell letter for l = 0..max_angular_momentum.
char angular_momentum_letter(int l);

// One contracted shell of an element as read from a basis set file:
// chi(r) = sum_k coefficients[k] * r^l exp(-exponents[k] r^2), not normalised.
// Coefficients are stored exactly as given in the file.
struct ShellDefinition {
  int l = 0;
  std::vector<double> exponents;
  std::vector<double> coefficients;

  std::size_t n_primitives() const { return exponents.size(); }
};

// A basis set for a range of elements, read from an NWChem-format file as
// provided by the Basis Set Exchange.
//
// Conventions (chosen to match the reference code PySCF, so that matrices
// can be compared element by element):
//
//  * Segmented representation. A generally contracted block with n
//    coefficient columns becomes n shells, in column order. A combined block
//    (e.g. "SP") becomes one shell per letter.
//  * Primitives whose coefficient is exactly zero in a column are dropped
//    from that shell.
//  * The shells of an element are ordered by angular momentum; within one
//    angular momentum the file order is kept.
class BasisLibrary {
 public:
  // Parses the NWChem format: a block `BASIS "<label>" [SPHERICAL|CARTESIAN]
  // [PRINT]` ... `END` with shell headers `<element> <type>` followed by
  // lines `<exponent> <coefficient> [<coefficient> ...]`. Lines starting with
  // '#' are comments. Elements beyond Ar are skipped. ECP blocks are rejected.
  // Throws std::runtime_error with the line number on malformed input.
  static BasisLibrary parse_nwchem(std::istream& in, std::string name);

  // Loads `<basis directory>/<file name for basis>`; see basis_file_name().
  static BasisLibrary load(std::string_view basis_name);
  static BasisLibrary load_file(const std::filesystem::path& path, std::string name);

  const std::string& name() const { return name_; }

  // True if the file declares spherical functions (the default if a molecule
  // basis is built without an explicit choice).
  bool spherical() const { return spherical_; }

  bool has_element(int Z) const { return shells_.count(Z) != 0; }

  // Shells of element Z. Throws std::out_of_range if the element is missing.
  const std::vector<ShellDefinition>& shells(int Z) const;

 private:
  std::string name_;
  bool spherical_ = false;
  std::map<int, std::vector<ShellDefinition>> shells_;
};

// File name of a basis set: lower case, '*' replaced by "_st", extension
// ".nw" (e.g. "6-31G*" -> "6-31g_st.nw"). Same convention as
// scripts/fetch_basis.py.
std::string basis_file_name(std::string_view basis_name);

// Directory with the basis set files: environment variable ELLIDFT_BASIS_DIR
// if set, otherwise the data/basis directory of the source tree.
std::filesystem::path basis_directory();

}  // namespace ellidft
