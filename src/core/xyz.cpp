#include "core/xyz.hpp"

#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>

#include "core/elements.hpp"
#include "core/parse_util.hpp"
#include "core/units.hpp"

namespace ellidft {

namespace {

[[noreturn]] void fail(int line, const std::string& what) {
  throw std::runtime_error("XYZ line " + std::to_string(line) + ": " + what);
}

}  // namespace

Molecule read_xyz(std::istream& in, LengthUnit unit) {
  std::string line;
  int line_no = 0;

  if (!std::getline(in, line)) fail(1, "missing atom count");
  ++line_no;
  const auto count_tokens = detail::split_whitespace(line);
  if (count_tokens.size() != 1) fail(line_no, "expected a single atom count");
  const auto count = detail::parse_unsigned(count_tokens[0]);
  if (!count || *count < 1) fail(line_no, "invalid atom count '" + count_tokens[0] + "'");

  if (!std::getline(in, line)) fail(2, "missing comment line");
  ++line_no;

  Molecule mol;
  mol.name = std::string(detail::trim(line));
  mol.atoms.reserve(static_cast<std::size_t>(*count));

  // Division (not multiplication by the reciprocal), so that the result is
  // bit-identical to `x / 0.529177210903` in the reference scripts.
  const double divisor = (unit == LengthUnit::Angstrom) ? units::bohr_in_angstrom : 1.0;

  for (long i = 0; i < *count; ++i) {
    if (!std::getline(in, line))
      fail(line_no + 1, "file ends after " + std::to_string(i) + " of " +
                            std::to_string(*count) + " atoms");
    ++line_no;
    const auto tok = detail::split_whitespace(line);
    if (tok.size() < 4) fail(line_no, "expected '<symbol> <x> <y> <z>'");

    Atom atom;
    atom.Z = atomic_number(tok[0]);
    if (atom.Z == 0)
      fail(line_no, "unknown or unsupported element '" + tok[0] + "' (supported: H to Ar)");
    for (std::size_t k = 0; k < 3; ++k) {
      const auto v = detail::parse_double(tok[k + 1]);
      if (!v) fail(line_no, "invalid coordinate '" + tok[k + 1] + "'");
      atom.r[k] = *v / divisor;
    }
    mol.atoms.push_back(atom);
  }

  // Only blank lines may follow; anything else means the count is wrong.
  while (std::getline(in, line)) {
    ++line_no;
    if (!detail::trim(line).empty())
      fail(line_no, "unexpected content after " + std::to_string(*count) + " atoms");
  }
  return mol;
}

Molecule read_xyz_file(const std::filesystem::path& path, LengthUnit unit) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open XYZ file '" + path.string() + "'");
  try {
    return read_xyz(in, unit);
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(path.string() + ": " + e.what());
  }
}

}  // namespace ellidft
