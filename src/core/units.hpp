#pragma once

// Physical constants and unit conversions.
//
// elliDFT uses atomic units (Hartree, Bohr) internally. Conversions happen
// only at input and output.
//
// Values: CODATA 2018, as fixed in the project specification
// ("Agenten-Anweisungen", section 5 "Einheiten und Konstanten" and section 8
// "Umrechnungen").

namespace ellidft::units {

// 1 Bohr in Angstrom (spec section 5).
inline constexpr double bohr_in_angstrom = 0.529177210903;

// 1 Hartree in eV (spec section 5).
inline constexpr double hartree_in_ev = 27.211386245988;

inline constexpr double angstrom_to_bohr(double x) { return x / bohr_in_angstrom; }
inline constexpr double bohr_to_angstrom(double x) { return x * bohr_in_angstrom; }

}  // namespace ellidft::units
