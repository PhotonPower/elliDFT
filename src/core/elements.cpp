#include "core/elements.hpp"

#include <array>
#include <cctype>
#include <stdexcept>

namespace ellidft {

namespace {

const std::array<std::string, max_atomic_number> symbols = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
    "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar"};

bool equal_ignore_case(std::string_view a, std::string_view b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i) {
    if (std::tolower(static_cast<unsigned char>(a[i])) !=
        std::tolower(static_cast<unsigned char>(b[i])))
      return false;
  }
  return true;
}

}  // namespace

int atomic_number(std::string_view symbol) {
  for (int z = 1; z <= max_atomic_number; ++z) {
    if (equal_ignore_case(symbol, symbols[static_cast<std::size_t>(z - 1)])) return z;
  }
  return 0;
}

const std::string& element_symbol(int Z) {
  if (Z < 1 || Z > max_atomic_number)
    throw std::out_of_range("element_symbol: unsupported atomic number " + std::to_string(Z));
  return symbols[static_cast<std::size_t>(Z - 1)];
}

}  // namespace ellidft
