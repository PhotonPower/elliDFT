#pragma once

#include <string>
#include <string_view>

namespace ellidft {

// Highest supported atomic number.
//
// The project scope is H to Cl (Z = 1..17). Argon (Z = 18) is included because
// the validation table of the specification (section 7) lists Ar as a test
// atom. Nothing beyond Ar is supported.
inline constexpr int max_atomic_number = 18;

// Atomic number for an element symbol ("H", "he", "CL", ...), matched without
// regard to case. Returns 0 if the symbol is unknown or beyond Ar.
int atomic_number(std::string_view symbol);

// Element symbol in standard capitalisation ("H", "He", "Cl").
// Throws std::out_of_range for Z outside 1..max_atomic_number.
const std::string& element_symbol(int Z);

}  // namespace ellidft
