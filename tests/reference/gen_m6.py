#!/usr/bin/env python3
"""Guide values for milestone M6: measured refractive indices of polystyrene and PMMA.

Writes data/m6/guide_values.json: n(lambda) from 400 to 1000 nm and the
polarisability volume per repeat unit that the Lorentz-Lorenz relation
assigns to n at 589 nm.

Sources
  Refractive index: Sellmeier fits in the refractiveindex.info database
  (M. N. Polyanskiy, Sci. Data 11, 94 (2024)), pages "PS" and "PMMA",
  consulted 2026-11-14. The database attributes the data to N. Sultanova,
  S. Kasarova, I. Nikolov, "Dispersion properties of optical polymers",
  Acta Physica Polonica A 116, 585-587 (2009); range 0.4368-1.052 um, 20 C.
  The original paper was not consulted.
      n^2 - 1 = B lambda^2 / (lambda^2 - C),  lambda in micrometres
      polystyrene: B = 1.4435, C = 0.020216
      PMMA:        B = 1.1819, C = 0.011313

  Density: 1.05 g/cm^3 (polystyrene) and 1.18 g/cm^3 (PMMA) are the rounded
  textbook values of the specification. NO SOURCE has been checked for them
  yet; the specification asks for one before M6 is accepted.

Does not use PySCF.

Usage: python tests/reference/gen_m6.py
"""

import json
import math
import pathlib

AVOGADRO_FACTOR = 0.602214076  # N [1/A^3] = factor * rho [g/cm^3] / M [g/mol]  (specification, section 8)

MATERIALS = {
    "polystyrene": {"repeat_unit": "C8H8", "molar_mass_g_mol": 104.15, "density_g_cm3": 1.05,
                    "sellmeier_B": 1.4435, "sellmeier_C_um2": 0.020216},
    "pmma": {"repeat_unit": "C5H8O2", "molar_mass_g_mol": 100.12, "density_g_cm3": 1.18,
             "sellmeier_B": 1.1819, "sellmeier_C_um2": 0.011313},
}
WAVELENGTHS_NM = [400, 450, 500, 550, 589, 600, 650, 700, 750, 800, 850, 900, 950, 1000]
# The fits cover 436.8 to 1052 nm; 400 nm is an extrapolation and marked as such.
FIT_RANGE_NM = (436.8, 1052.0)


def refractive_index(material: dict, wavelength_nm: float) -> float:
    lam2 = (wavelength_nm / 1000.0) ** 2
    return math.sqrt(1.0 + material["sellmeier_B"] * lam2 / (lam2 - material["sellmeier_C_um2"]))


def polarizability_volume(n: float, density: float, molar_mass: float) -> float:
    """Lorentz-Lorenz: (n^2 - 1)/(n^2 + 2) = (4 pi / 3) N alpha'."""
    number_density = AVOGADRO_FACTOR * density / molar_mass
    return (n * n - 1.0) / (n * n + 2.0) * 3.0 / (4.0 * math.pi * number_density)


def main() -> None:
    out = {
        "meta": {
            "generated_by": "tests/reference/gen_m6.py",
            "refractive_index_source": "refractiveindex.info (Polyanskiy, Sci. Data 11, 94 (2024)), Sellmeier fits "
                                       "to Sultanova, Kasarova, Nikolov, Acta Phys. Pol. A 116, 585 (2009); 20 C; "
                                       "database pages consulted 2026-11-14, original paper not consulted",
            "density_source": "NONE CHECKED: rounded textbook values taken from the specification",
            "fit_range_nm": list(FIT_RANGE_NM),
        },
        "materials": {},
    }
    for name, material in MATERIALS.items():
        n589 = refractive_index(material, 589.0)
        out["materials"][name] = {
            **material,
            "n_589nm": n589,
            "polarizability_volume_per_unit_A3": polarizability_volume(
                n589, material["density_g_cm3"], material["molar_mass_g_mol"]),
            "refractive_index": [
                {"wavelength_nm": w, "n": refractive_index(material, w),
                 "extrapolated": not (FIT_RANGE_NM[0] <= w <= FIT_RANGE_NM[1])}
                for w in WAVELENGTHS_NM],
        }
        print(f"{name}: n(589 nm) = {n589:.5f}, alpha' = {out['materials'][name]['polarizability_volume_per_unit_A3']:.3f} A^3")
    path = pathlib.Path(__file__).resolve().parent / "data" / "m6" / "guide_values.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {path.name}")


if __name__ == "__main__":
    main()
