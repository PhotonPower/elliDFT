"""Wall time of a PySCF ground-state calculation, as the reference for the
run-time criteria of the specification (M3: at most ten times the reference).

    python tests/reference/time_scf.py <molecule> <basis> <functional> [--aux <basis>] [--grid 75,302]

<molecule> is a name from data/geometries. Same geometry files, basis set
files, unpruned grid and thresholds as the other reference scripts. Prints
one JSON object.
"""
import argparse
import json
import sys
import time

import pyscf
from pyscf import dft, lib, scf

import common
import gen_m2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("molecule")
    parser.add_argument("basis")
    parser.add_argument("functional")
    parser.add_argument("--aux", default=None, help="auxiliary basis for density fitting")
    parser.add_argument("--grid", default="75,302")
    args = parser.parse_args()
    n_radial, n_angular = (int(x) for x in args.grid.split(","))

    start = time.perf_counter()
    symbols, coords, _ = common.read_xyz_bohr(args.molecule)
    mol = common.build_mol(symbols, coords, args.basis, cart=not common.basis_is_spherical(args.basis))
    if args.functional.upper() == "HF":
        mf = scf.RHF(mol)
    else:
        mf = dft.RKS(mol)
        mf.xc = args.functional
        mf.grids = gen_m2.make_grid(mol, n_radial, n_angular)
    if args.aux:
        mf = mf.density_fit(auxbasis=common.pyscf_basis(args.aux, symbols))
    mf.conv_tol = 1e-10
    mf.verbose = 0
    energy = mf.kernel()
    seconds = time.perf_counter() - start
    json.dump({"program": "PySCF " + pyscf.__version__, "molecule": args.molecule, "basis": args.basis,
               "functional": args.functional, "aux_basis": args.aux, "grid": [n_radial, n_angular],
               "threads": lib.num_threads(), "converged": bool(mf.converged), "energy_hartree": float(energy),
               "iterations": int(getattr(mf, "cycles", -1)), "wall_seconds": seconds}, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
