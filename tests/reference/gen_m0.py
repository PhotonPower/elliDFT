#!/usr/bin/env python3
"""Reference data for milestone M0: basis set structure as seen by PySCF.

Writes
  data/m0/basis_shells.json  shells of every element H..Ar in every basis set
  data/m0/molecules.json     geometry in Bohr and basis sizes of all test molecules

Usage: python tests/reference/gen_m0.py
"""

import common


def element_shells() -> dict:
    out = {}
    for basis_name in common.ALL_BASIS_SETS:
        elements = {}
        for symbol in common.ELEMENTS:
            # A single atom; the electronic state is irrelevant for the basis.
            mol = common.build_mol([symbol], [[0.0, 0.0, 0.0]], basis_name,
                                   cart=not common.basis_is_spherical(basis_name))
            elements[symbol] = [{k: s[k] for k in ("l", "exponents", "coefficients")}
                                for s in common.expanded_shells(mol)]
        out[basis_name] = {
            "file": common.basis_file(basis_name).name,
            "spherical": common.basis_is_spherical(basis_name),
            "elements": elements,
        }
    return out


def basis_size(symbols, coords, basis_name: str) -> dict:
    sizes = {}
    for label, cart in (("spherical", False), ("cartesian", True)):
        mol = common.build_mol(symbols, coords, basis_name, cart=cart)
        shells = common.expanded_shells(mol)
        per_l = [0] * (max(s["l"] for s in shells) + 1)
        for s in shells:
            per_l[s["l"]] += 1
        assert len(shells) == int(sum(mol.bas_nctr(i) for i in range(mol.nbas)))
        sizes[label] = {
            "n_functions": int(mol.nao),            # PySCF: mol.nao
            "n_shells": len(shells),                # contracted shells, segmented count
            "n_shells_per_l": per_l,
            "n_primitives": sum(len(s["exponents"]) for s in shells),
            "pyscf_nbas": int(mol.nbas),            # general contractions counted once
        }
    return sizes


def molecules() -> dict:
    out = {}
    for stem in common.ALL_MOLECULES:
        symbols, coords, comment = common.read_xyz_bohr(stem)
        probe = common.build_mol(symbols, coords, "STO-3G", cart=False)
        out[stem] = {
            "xyz_file": f"data/geometries/{stem}.xyz",
            "comment": comment,
            "charge": 0,
            "n_electrons": int(probe.nelectron),
            "symbols": symbols,
            "geometry_bohr": coords,
            "basis": {name: basis_size(symbols, coords, name) for name in common.ALL_BASIS_SETS},
        }
        print(f"  {stem}: {len(symbols)} atoms")
    return out


def main() -> None:
    info = common.meta("gen_m0.py", basis_source="data/basis/*.nw (Basis Set Exchange)",
                       basis_parser="pyscf.gto.basis.parse_nwchem.load(optimize=False)")
    common.write_json(common.DATA_DIR / "m0" / "basis_shells.json",
                      {"meta": info, "basis_sets": element_shells()})
    common.write_json(common.DATA_DIR / "m0" / "molecules.json",
                      {"meta": info, "molecules": molecules()})


if __name__ == "__main__":
    main()
