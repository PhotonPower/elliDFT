#!/usr/bin/env python3
"""Reference data for milestone M4: hybrid functionals and erf-attenuated integrals (PySCF).

Writes to data/m4/:
  erf_ints_<molecule>_<basis>.json  two-electron integrals with erf(omega r)/r:
                                    four-centre (8-fold packed), three- and two-centre
  hybrids.json                      PBE0, B3LYP and CAM-B3LYP energies on the (99, 590) grid

Functionals are given by their libxc names, so that the B3LYP variant is
unambiguous: HYB_GGA_XC_B3LYP uses VWN-RPA.

Usage: python tests/reference/gen_m4.py [--only ints|energies]
"""

import argparse

import numpy
from pyscf import df, dft

import common
from gen_m2 import make_grid

OUT = common.DATA_DIR / "m4"
OMEGA = 0.33
AUX_BASIS = "def2-universal-JKFIT"
KS_GRID = (99, 590)
FUNCTIONALS = {"PBE0": "HYB_GGA_XC_PBEH", "B3LYP": "HYB_GGA_XC_B3LYP", "CAM-B3LYP": "HYB_GGA_XC_CAM_B3LYP"}
ENERGY_CASES = [(stem, "def2-SVP") for stem in common.SMALL_MOLECULES]


def gen_ints() -> None:
    for stem, basis_name in (("h2o", "STO-3G"), ("h2o", "def2-SVP"), ("hcl", "STO-3G")):
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        aux = common.build_mol(symbols, coords, AUX_BASIS, cart=False)
        n = mol.nao
        with mol.with_range_coulomb(OMEGA):
            four = mol.intor("int2e", aosym="s8")
        coulomb = mol.intor("int2e", aosym="s8")
        assert numpy.abs(four).max() < numpy.abs(coulomb).max()  # attenuated integrals are smaller
        with aux.with_range_coulomb(OMEGA):
            two = aux.intor("int2c2e")
        with mol.with_range_coulomb(OMEGA), aux.with_range_coulomb(OMEGA):
            three = df.incore.aux_e2(mol, aux, intor="int3c2e", aosym="s1")
        three_coulomb = df.incore.aux_e2(mol, aux, intor="int3c2e", aosym="s1")
        assert numpy.abs(three - three_coulomb).max() > 1e-3  # the range separation took effect
        rows, cols = numpy.tril_indices(n)
        common.write_json(OUT / f"erf_ints_{stem}_{basis_name.lower()}.json", {
            "meta": common.meta("gen_m4.py", operator="erf(omega r12)/r12", omega=OMEGA,
                                auxiliary_basis=AUX_BASIS),
            "molecule": {"charge": 0, "symbols": symbols, "geometry_bohr": coords},
            "basis": {"name": basis_name, "spherical": spherical},
            "auxiliary_basis": AUX_BASIS,
            "omega": OMEGA,
            "n_functions": int(n),
            "n_auxiliary": int(aux.nao),
            "eri_packed": four.tolist(),
            "two_center": two[numpy.tril_indices(aux.nao)].tolist(),
            "three_center": three[rows, cols, :].tolist() if n <= 12 else None,
        })


def gen_energies() -> None:
    settings = {
        "functionals": FUNCTIONALS,
        "grid": {"n_radial": KS_GRID[0], "n_angular": KS_GRID[1]},
        "small_rho_cutoff": 0.0,
        "conv_tol": 1e-12,
        "conv_tol_grad": 1e-9,
        "exchange": "exact four-centre integrals (no density fitting)",
    }
    records = []
    coefficients = {}
    for stem, basis_name in ENERGY_CASES:
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        for label, xc in FUNCTIONALS.items():
            mf = dft.RKS(mol)
            mf.xc = xc
            mf.grids = make_grid(mol, *KS_GRID)
            mf.small_rho_cutoff = 0.0
            mf.conv_tol = settings["conv_tol"]
            mf.conv_tol_grad = settings["conv_tol_grad"]
            mf.direct_scf_tol = 1e-14
            mf.max_cycle = 200
            mf.kernel()
            assert mf.converged, (stem, label)
            omega, alpha, hyb = mf._numint.rsh_and_hybrid_coeff(mf.xc)
            # PySCF: hyb = short-range fraction, alpha = long-range fraction of exact exchange.
            coefficients[label] = {"omega": float(omega), "exact_exchange": float(hyb),
                                   "long_range_exchange": float(alpha - hyb) if omega != 0 else 0.0}
            records.append({
                "molecule": stem, "basis": basis_name, "spherical": spherical, "functional": label,
                "charge": 0, "symbols": symbols, "geometry_bohr": coords,
                "n_electrons": int(mol.nelectron), "n_functions": int(mol.nao),
                "energy_total": float(mf.e_tot),
                "dipole_au": [float(x) for x in mf.dip_moment(unit="AU", verbose=0)],
                "scf_cycles": int(mf.cycles),
            })
            print(f"  {stem:8s} {label:10s} E = {mf.e_tot:.10f}  cycles {mf.cycles}", flush=True)
    settings["exchange_coefficients_in_pyscf"] = coefficients
    common.write_json(OUT / "hybrids.json", {"meta": common.meta("gen_m4.py", **settings), "results": records})


def main() -> None:
    parts = {"ints": gen_ints, "energies": gen_energies}
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", choices=sorted(parts), help="generate one part only")
    args = parser.parse_args()
    for name, function in parts.items():
        if args.only in (None, name):
            function()


if __name__ == "__main__":
    main()
