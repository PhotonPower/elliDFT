#pragma once

// Small text parsing helpers shared by the input readers.

#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace ellidft::detail {

inline std::string_view trim(std::string_view s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.remove_prefix(1);
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.remove_suffix(1);
  return s;
}

inline std::vector<std::string> split_whitespace(std::string_view s) {
  std::vector<std::string> out;
  std::size_t i = 0;
  while (i < s.size()) {
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    std::size_t j = i;
    while (j < s.size() && !std::isspace(static_cast<unsigned char>(s[j]))) ++j;
    if (j > i) out.emplace_back(s.substr(i, j - i));
    i = j;
  }
  return out;
}

// Parses a floating point number. The whole token must be consumed.
// With `fortran_exponent`, 'D'/'d' is accepted as exponent marker (1.0D-03).
inline std::optional<double> parse_double(std::string token, bool fortran_exponent = false) {
  if (token.empty()) return std::nullopt;
  if (fortran_exponent) {
    for (char& c : token)
      if (c == 'D' || c == 'd') c = 'e';
  }
  // Reject things strtod accepts but that are not plain decimal numbers
  // (hex floats, "inf", "nan").
  for (char c : token) {
    const bool ok = std::isdigit(static_cast<unsigned char>(c)) || c == '+' || c == '-' ||
                    c == '.' || c == 'e' || c == 'E';
    if (!ok) return std::nullopt;
  }
  errno = 0;
  char* end = nullptr;
  const double v = std::strtod(token.c_str(), &end);
  if (end != token.c_str() + token.size() || errno == ERANGE) return std::nullopt;
  return v;
}

// Parses a non-negative decimal integer. The whole token must be consumed.
inline std::optional<long> parse_unsigned(const std::string& token) {
  if (token.empty()) return std::nullopt;
  for (char c : token)
    if (!std::isdigit(static_cast<unsigned char>(c))) return std::nullopt;
  errno = 0;
  const long v = std::strtol(token.c_str(), nullptr, 10);
  if (errno == ERANGE) return std::nullopt;
  return v;
}

}  // namespace ellidft::detail
