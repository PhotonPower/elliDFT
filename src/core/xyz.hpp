#pragma once

#include <filesystem>
#include <iosfwd>

#include "core/molecule.hpp"

namespace ellidft {

enum class LengthUnit { Angstrom, Bohr };

// Reads a molecule in XYZ format:
//
//   line 1:  number of atoms N
//   line 2:  comment (stored as Molecule::name, surrounding blanks removed)
//   N lines: <element symbol> <x> <y> <z>
//
// Coordinates are interpreted in `unit` (XYZ files are in Angstrom by
// convention) and stored in Bohr. Extra columns after z are ignored; blank
// lines after the last atom are allowed. Throws std::runtime_error with the
// line number on malformed input or unsupported elements.
Molecule read_xyz(std::istream& in, LengthUnit unit = LengthUnit::Angstrom);

Molecule read_xyz_file(const std::filesystem::path& path,
                       LengthUnit unit = LengthUnit::Angstrom);

}  // namespace ellidft
