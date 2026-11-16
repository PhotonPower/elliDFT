#!/usr/bin/env python3
"""Reference data for milestone M3: density fitting (PySCF).

Writes to data/m3/:
  df_ints_<molecule>.json   two- and three-centre Coulomb integrals
  df.json                   Hartree-Fock and PBE energies with density fitting,
                            together with the energies without it

Density fitting in the reference: Coulomb metric, auxiliary basis read from
the same file as in elliDFT (def2-universal-JKFIT), spherical functions. The
energy cases use spherical orbital basis sets, for which PySCF's own choice of
the auxiliary function type is spherical as well (asserted below).

Usage: python tests/reference/gen_m3.py [--only ints|energies]
"""

import argparse

import numpy
from pyscf import df, dft, scf

import common
from gen_m2 import make_grid

OUT = common.DATA_DIR / "m3"
AUX_BASIS = "def2-universal-JKFIT"
KS_GRID = (99, 590)
ENERGY_CASES = [(stem, "def2-SVP") for stem in common.SMALL_MOLECULES]


def aux_mol(symbols, coords):
    """Auxiliary basis as its own molecule, always with spherical functions.

    PySCF's make_auxmol would copy the function type of the orbital basis
    (cartesian for 6-31G*); elliDFT uses the auxiliary basis as defined."""
    return common.build_mol(symbols, coords, AUX_BASIS, cart=False)


def gen_ints() -> None:
    # Spherical orbital basis sets only: PySCF cannot combine a cartesian
    # orbital basis (6-31G*) with a spherical auxiliary basis.
    for stem, basis_name in (("h2o", "def2-SVP"), ("hcl", "def2-SVP")):
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        aux = aux_mol(symbols, coords)
        two = aux.intor("int2c2e")
        three = df.incore.aux_e2(mol, aux, intor="int3c2e", aosym="s1")  # (nao, nao, naux)
        n = mol.nao
        rows, cols = numpy.tril_indices(n)
        common.write_json(OUT / f"df_ints_{stem}.json", {
            "meta": common.meta("gen_m3.py", auxiliary_basis=AUX_BASIS, metric="Coulomb",
                                packing_three_center="rows (mu nu) with mu >= nu, index mu(mu+1)/2 + nu; columns P"),
            "molecule": {"charge": 0, "symbols": symbols, "geometry_bohr": coords},
            "basis": {"name": basis_name, "spherical": spherical},
            "auxiliary_basis": AUX_BASIS,
            "n_functions": int(n),
            "n_auxiliary": int(aux.nao),
            "two_center": two[numpy.tril_indices(aux.nao)].tolist(),
            "three_center": three[rows, cols, :].tolist(),
        })


def gen_energies() -> None:
    settings = {
        "auxiliary_basis": AUX_BASIS,
        "metric": "Coulomb",
        "conv_tol": 1e-12,
        "conv_tol_grad": 1e-9,
        "grid": {"n_radial": KS_GRID[0], "n_angular": KS_GRID[1]},
        "small_rho_cutoff": 0.0,
        "functionals": {"HF": "exact exchange with RI-J and RI-K", "PBE": "GGA_X_PBE,GGA_C_PBE with RI-J"},
    }
    records = []
    for stem, basis_name in ENERGY_CASES:
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        auxbasis = common.pyscf_basis(AUX_BASIS, symbols)
        for label in ("HF", "PBE"):
            energies = {}
            for fitted in (True, False):
                if label == "HF":
                    mf = scf.RHF(mol)
                else:
                    mf = dft.RKS(mol)
                    mf.xc = "GGA_X_PBE,GGA_C_PBE"
                    mf.grids = make_grid(mol, *KS_GRID)
                    mf.small_rho_cutoff = 0.0
                if fitted:
                    mf = mf.density_fit(auxbasis=auxbasis)
                mf.conv_tol = settings["conv_tol"]
                mf.conv_tol_grad = settings["conv_tol_grad"]
                mf.direct_scf_tol = 1e-14
                mf.max_cycle = 200
                mf.kernel()
                assert mf.converged, (stem, label, fitted)
                energies[fitted] = float(mf.e_tot)
                if fitted:
                    n_aux = int(mf.with_df.get_naoaux())
                    assert not mf.with_df.auxmol.cart
            records.append({
                "molecule": stem, "basis": basis_name, "spherical": spherical, "functional": label,
                "charge": 0, "symbols": symbols, "geometry_bohr": coords,
                "n_electrons": int(mol.nelectron), "n_functions": int(mol.nao), "n_auxiliary": n_aux,
                "energy_total": energies[True],
                "energy_total_without_fitting": energies[False],
                "fitting_error": energies[True] - energies[False],
            })
            print(f"  {stem:8s} {label:4s} E(DF) = {energies[True]:.10f}  DF error = "
                  f"{energies[True] - energies[False]:+.2e}", flush=True)
    common.write_json(OUT / "df.json", {"meta": common.meta("gen_m3.py", **settings), "results": records})


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
