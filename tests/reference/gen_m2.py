#!/usr/bin/env python3
"""Reference data for milestone M2: grids, functionals, Kohn-Sham energies (PySCF).

Writes to data/m2/:
  grid_<molecule>.json      small molecular grids, point by point
  ks.json                   LDA and PBE total energies on the (99, 590) grid

Grid settings of the reference, set explicitly (they are also PySCF's
defaults except for pruning):
  radial grid       Treutler-Ahlrichs (M4), element-specific radius parameter
  angular grid      Lebedev, the same number of points on every radial shell
  pruning           none
  partitioning      Becke, with the size adjustment of Treutler and Ahlrichs
                    based on Bragg radii

Usage: python tests/reference/gen_m2.py [--only grid|ks]
"""

import argparse

import numpy
from pyscf import dft
from pyscf.dft import gen_grid, radi

import common

OUT = common.DATA_DIR / "m2"

GRID_SETTINGS = {
    "radial": "Treutler-Ahlrichs M4 (pyscf.dft.radi.treutler), ATOM_SPECIFIC_TREUTLER_GRIDS = True",
    "angular": "Lebedev, unpruned",
    "partition": "Becke (original_becke) with treutler_atomic_radii_adjust, Bragg radii",
}


def make_grid(mol, n_radial: int, n_angular: int) -> gen_grid.Grids:
    assert radi.ATOM_SPECIFIC_TREUTLER_GRIDS
    grids = gen_grid.Grids(mol)
    grids.atom_grid = (n_radial, n_angular)
    grids.radi_method = radi.treutler
    grids.atomic_radii = radi.BRAGG_RADII
    grids.radii_adjust = radi.treutler_atomic_radii_adjust
    grids.becke_scheme = gen_grid.original_becke
    grids.prune = None
    return grids


def gen_grid_cases() -> None:
    # Small grids, stored point by point. H2O: two elements with different
    # radii; HCl: third period; CH4: five atoms.
    for stem, n_radial, n_angular in (("h2o", 12, 26), ("hcl", 10, 14), ("ch4", 8, 14)):
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, "STO-3G", cart=False)
        grids = make_grid(mol, n_radial, n_angular)
        grids.build(with_non0tab=False, sort_grids=False)
        n = mol.natm * n_radial * n_angular  # PySCF appends zero-weight padding points
        assert numpy.all(grids.weights[n:] == 0.0)
        common.write_json(OUT / f"grid_{stem}.json", {
            "meta": common.meta("gen_m2.py", **GRID_SETTINGS),
            "molecule": {"charge": 0, "symbols": symbols, "geometry_bohr": coords},
            "n_radial": n_radial,
            "n_angular": n_angular,
            "points": grids.coords[:n].tolist(),
            "weights": grids.weights[:n].tolist(),
            "atom": grids.atm_idx[:n].tolist(),
        })


def gen_basis_values() -> None:
    """Basis functions and gradients at a few points (pyscf eval_gto)."""
    points = [[0.0, 0.0, 0.0], [0.3, -0.2, 0.5], [-1.1, 0.7, 0.4], [0.05, 1.4, 1.1],
              [2.5, -1.9, 0.8], [0.0, 0.0, 6.0], [-0.4, -0.6, -2.2], [1.0e-3, 2.0e-3, -1.0e-3]]
    cases = []
    for stem, basis_name, spherical in (("h2o", "def2-TZVP", True), ("h2o", "def2-TZVP", False),
                                        ("hcl", "6-31G*", False), ("h2o", "def2-universal-JKFIT", True)):
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        name = "GTOval_sph_deriv1" if spherical else "GTOval_cart_deriv1"
        values = mol.eval_gto(name, numpy.array(points))  # (4, n_points, n_functions): value, d/dx, d/dy, d/dz
        cases.append({
            "molecule": {"charge": 0, "symbols": symbols, "geometry_bohr": coords},
            "basis": {"name": basis_name, "spherical": spherical},
            "n_functions": int(mol.nao),
            "value": values[0].tolist(), "dx": values[1].tolist(),
            "dy": values[2].tolist(), "dz": values[3].tolist(),
        })
    common.write_json(OUT / "basis_values.json", {
        "meta": common.meta("gen_m2.py", evaluation="mol.eval_gto GTOval_{sph,cart}_deriv1"),
        "points_bohr": points,
        "cases": cases,
    })


# Functionals by explicit libxc names, so that the VWN variant is unambiguous
# (LDA_C_VWN is VWN5).
KS_FUNCTIONALS = {"LDA": "LDA_X,LDA_C_VWN", "PBE": "GGA_X_PBE,GGA_C_PBE"}
KS_GRID = (99, 590)
KS_CASES = [(stem, "def2-SVP") for stem in common.SMALL_MOLECULES] + \
           [("h2o", "6-31G*"), ("hcl", "6-31G*"), ("h2o", "cc-pVDZ"), ("ne", "aug-cc-pVDZ")]


def gen_ks() -> None:
    settings = {
        "method": "RKS, closed shell, exact four-centre Coulomb integrals (no density fitting)",
        "functionals": KS_FUNCTIONALS,
        "grid": {"n_radial": KS_GRID[0], "n_angular": KS_GRID[1], **GRID_SETTINGS},
        "small_rho_cutoff": 0.0,
        "conv_tol": 1e-12,
        "conv_tol_grad": 1e-9,
        "direct_scf_tol": 1e-14,
        "init_guess": "minao",
        "libxc_version_of_pyscf": dft.libxc.__version__,
    }
    records = []
    for stem, basis_name in KS_CASES:
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        for label, xc in KS_FUNCTIONALS.items():
            mf = dft.RKS(mol)
            mf.xc = xc
            mf.grids = make_grid(mol, *KS_GRID)
            mf.small_rho_cutoff = settings["small_rho_cutoff"]  # keep every grid point
            mf.conv_tol = settings["conv_tol"]
            mf.conv_tol_grad = settings["conv_tol_grad"]
            mf.direct_scf_tol = settings["direct_scf_tol"]
            mf.init_guess = settings["init_guess"]
            mf.max_cycle = 200
            mf.kernel()
            assert mf.converged, (stem, basis_name, label)
            assert not hasattr(mf, "with_df"), "density fitting must be off"
            dm = mf.make_rdm1()
            n_grid, e_xc, _ = mf._numint.nr_rks(mol, mf.grids, mf.xc, dm)
            records.append({
                "molecule": stem,
                "basis": basis_name,
                "spherical": spherical,
                "functional": label,
                "charge": 0,
                "symbols": symbols,
                "geometry_bohr": coords,
                "n_electrons": int(mol.nelectron),
                "n_functions": int(mol.nao),
                "n_grid_points": int((mf.grids.weights != 0).sum()),
                "energy_total": float(mf.e_tot),
                "energy_xc": float(e_xc),
                "energy_nuclear_repulsion": float(mf.energy_nuc()),
                "electrons_on_grid": float(n_grid),
                "dipole_au": [float(x) for x in mf.dip_moment(unit="AU", verbose=0)],
                "homo": float(mf.mo_energy[mol.nelectron // 2 - 1]),
                "scf_cycles": int(mf.cycles),
            })
            print(f"  {stem:8s} {basis_name:12s} {label:4s} E = {mf.e_tot:.10f}  N(grid) - N = "
                  f"{n_grid - mol.nelectron:+.2e}  cycles {mf.cycles}", flush=True)
    common.write_json(OUT / "ks.json", {"meta": common.meta("gen_m2.py", **settings), "results": records})


def main() -> None:
    parts = {"grid": gen_grid_cases, "basis": gen_basis_values, "ks": gen_ks}
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", choices=sorted(parts), help="generate one part only")
    args = parser.parse_args()
    for name, function in parts.items():
        if args.only in (None, name):
            function()


if __name__ == "__main__":
    main()
