#include "basis/basis_library.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <stdexcept>

#include "core/elements.hpp"
#include "core/parse_util.hpp"

namespace ellidft {

namespace {

constexpr std::string_view shell_letters = "spdfghi";

std::string to_upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

bool is_alphabetic(const std::string& s) {
  return !s.empty() && std::all_of(s.begin(), s.end(), [](char c) {
    return std::isalpha(static_cast<unsigned char>(c)) != 0;
  });
}

[[noreturn]] void fail(int line, const std::string& what) {
  throw std::runtime_error("basis file line " + std::to_string(line) + ": " + what);
}

// A shell block while its data lines are being collected.
struct OpenBlock {
  int Z = 0;                 // 0: element beyond Ar, block is skipped
  std::vector<int> ls;       // one entry for a combined block per letter,
                             // a single entry otherwise
  bool combined = false;
  int header_line = 0;
  std::vector<std::vector<double>> rows;  // exponent followed by coefficients
};

void close_block(const OpenBlock& block, std::map<int, std::vector<ShellDefinition>>& out) {
  if (block.ls.empty()) return;
  if (block.rows.empty()) fail(block.header_line, "shell without primitives");
  if (block.Z == 0) return;

  const std::size_t n_columns = block.rows.front().size() - 1;
  if (block.combined && n_columns != block.ls.size())
    fail(block.header_line, "combined shell needs one coefficient column per letter");

  for (std::size_t col = 0; col < n_columns; ++col) {
    ShellDefinition shell;
    shell.l = block.combined ? block.ls[col] : block.ls.front();
    for (const auto& row : block.rows) {
      const double coefficient = row[col + 1];
      if (coefficient == 0.0) continue;  // primitive not part of this contraction
      shell.exponents.push_back(row[0]);
      shell.coefficients.push_back(coefficient);
    }
    if (shell.exponents.empty())
      fail(block.header_line, "contraction with only zero coefficients");
    out[block.Z].push_back(std::move(shell));
  }
}

}  // namespace

int angular_momentum_from_letter(char letter) {
  const char lower = static_cast<char>(std::tolower(static_cast<unsigned char>(letter)));
  const auto pos = shell_letters.find(lower);
  return pos == std::string_view::npos ? -1 : static_cast<int>(pos);
}

char angular_momentum_letter(int l) {
  if (l < 0 || l > max_angular_momentum)
    throw std::out_of_range("angular_momentum_letter: l = " + std::to_string(l));
  return shell_letters[static_cast<std::size_t>(l)];
}

BasisLibrary BasisLibrary::parse_nwchem(std::istream& in, std::string name) {
  BasisLibrary lib;
  lib.name_ = std::move(name);

  bool seen_basis = false;
  bool in_basis = false;
  OpenBlock block;
  std::string line;
  int line_no = 0;

  while (std::getline(in, line)) {
    ++line_no;
    const std::string_view text = detail::trim(line);
    if (text.empty() || text.front() == '#') continue;
    const auto tok = detail::split_whitespace(text);
    const std::string first = to_upper(tok[0]);

    if (!in_basis) {
      if (first == "BASIS") {
        if (seen_basis) fail(line_no, "more than one BASIS block");
        seen_basis = in_basis = true;
        // NWChem default without keyword is cartesian.
        lib.spherical_ = false;
        for (std::size_t i = 1; i < tok.size(); ++i) {
          const std::string key = to_upper(tok[i]);
          if (key == "SPHERICAL") lib.spherical_ = true;
          if (key == "CARTESIAN") lib.spherical_ = false;
        }
        continue;
      }
      if (first == "ECP")
        fail(line_no, "effective core potentials are not supported");
      fail(line_no, "unexpected content outside of a BASIS block");
    }

    if (first == "END") {
      close_block(block, lib.shells_);
      block = OpenBlock{};
      in_basis = false;
      continue;
    }

    if (is_alphabetic(tok[0])) {
      // Shell header: <element symbol> <shell type>
      close_block(block, lib.shells_);
      block = OpenBlock{};
      block.header_line = line_no;
      if (tok.size() != 2 || !is_alphabetic(tok[1]))
        fail(line_no, "expected '<element> <shell type>'");
      block.Z = atomic_number(tok[0]);  // 0 for elements beyond Ar: skipped
      block.combined = tok[1].size() > 1;
      for (char letter : tok[1]) {
        const int l = angular_momentum_from_letter(letter);
        if (l < 0) fail(line_no, "unknown shell type '" + tok[1] + "'");
        block.ls.push_back(l);
      }
      continue;
    }

    // Data line: exponent and one coefficient per contraction.
    if (block.ls.empty()) fail(line_no, "primitive data before any shell header");
    if (tok.size() < 2) fail(line_no, "expected exponent and at least one coefficient");
    std::vector<double> row;
    row.reserve(tok.size());
    for (const std::string& t : tok) {
      const auto v = detail::parse_double(t, /*fortran_exponent=*/true);
      if (!v) fail(line_no, "invalid number '" + t + "'");
      row.push_back(*v);
    }
    if (!(row[0] > 0.0)) fail(line_no, "exponent must be positive");
    if (!block.rows.empty() && block.rows.front().size() != row.size())
      fail(line_no, "inconsistent number of coefficient columns within a shell");
    block.rows.push_back(std::move(row));
  }

  if (in_basis) fail(line_no, "missing END of BASIS block");
  if (!seen_basis) fail(line_no, "no BASIS block found");

  // Order the shells of each element by angular momentum, keeping the file
  // order within one angular momentum (convention of the reference code).
  for (auto& [Z, shells] : lib.shells_) {
    std::stable_sort(shells.begin(), shells.end(),
                     [](const ShellDefinition& a, const ShellDefinition& b) { return a.l < b.l; });
  }
  return lib;
}

BasisLibrary BasisLibrary::load_file(const std::filesystem::path& path, std::string name) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open basis set file '" + path.string() + "'");
  try {
    return parse_nwchem(in, std::move(name));
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(path.string() + ": " + e.what());
  }
}

BasisLibrary BasisLibrary::load(std::string_view basis_name) {
  return load_file(basis_directory() / basis_file_name(basis_name), std::string(basis_name));
}

const std::vector<ShellDefinition>& BasisLibrary::shells(int Z) const {
  const auto it = shells_.find(Z);
  if (it == shells_.end())
    throw std::out_of_range("basis set '" + name_ + "' has no functions for Z = " +
                            std::to_string(Z));
  return it->second;
}

std::string basis_file_name(std::string_view basis_name) {
  std::string out;
  for (char c : basis_name) {
    if (c == '*')
      out += "_st";
    else
      out += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  }
  return out + ".nw";
}

std::filesystem::path basis_directory() {
  if (const char* env = std::getenv("ELLIDFT_BASIS_DIR"); env != nullptr && *env != '\0')
    return env;
  return ELLIDFT_DEFAULT_BASIS_DIR;
}

}  // namespace ellidft
