#!/usr/bin/env python3
"""Reference data for milestone M1: integrals and Hartree-Fock (PySCF).

Writes to data/m1/:
  ints1e_<molecule>_<basis>_<sph|cart>.json  overlap, kinetic, nuclear attraction, dipole
  eri_<molecule>_<basis>_<sph|cart>.json     two-electron integrals, 8-fold packed
  highl_<sph|cart>.json                      artificial basis with l = 0..6: one-electron
                                             matrices and selected two-electron shell quartets
  hf.json                                    restricted Hartree-Fock energies and dipoles

Conventions of the reference (PySCF / libcint), to be matched by elliDFT:

* Function order in a shell. Spherical: p as x, y, z; l >= 2 as m = -l..+l
  (real solid harmonics). Cartesian: lexicographic, xx xy xz yy yz zz, ...
* Normalisation. Spherical functions are normalised to 1. Cartesian functions
  with l >= 2 carry the normalisation of the radial part only, so the overlap
  diagonal is not 1 for all of them (see "overlap" in the cartesian files).
* nuclear = - sum_A Z_A <mu| 1/|r - R_A| |nu>, point nuclei.
* dipole_k = <mu| r_k |nu>, origin (0, 0, 0), no sign for the electron charge.
* Two-electron integrals (mu nu|lambda sigma) in chemists' notation.
* The dipole moment of the molecule is in atomic units (e*Bohr), nuclei minus
  electrons.

Usage: python tests/reference/gen_m1.py [--only ints1e|eri|highl|hf]
"""

import argparse

import numpy
from pyscf import scf
from pyscf.gto.basis import parse_nwchem

import common

OUT = common.DATA_DIR / "m1"

# Basis sets of the Hartree-Fock acceptance criterion (specification, M1).
HF_BASIS_SETS = ["STO-3G", "6-31G*", "def2-SVP", "cc-pVDZ", "aug-cc-pVDZ"]

# One-electron matrices are stored for all small molecules in these basis sets
# up to this size; larger cases are covered through the Hartree-Fock energies.
MAX_FUNCTIONS_1E = 60

# Additional one-electron cases: f functions in a regular basis set, both types.
EXTRA_1E = [("h2o", "def2-TZVP", False), ("h2o", "def2-TZVP", True),
            ("h2o", "6-31G*", False), ("hcl", "cc-pVDZ", True)]

# Full two-electron tensors (8-fold packed), "small molecules".
ERI_CASES = [("h2", "STO-3G"), ("h2", "def2-SVP"), ("he", "aug-cc-pVDZ"),
             ("h2o", "STO-3G"), ("h2o", "6-31G*"), ("h2o", "def2-SVP"),
             ("hf", "STO-3G"), ("nh3", "STO-3G"), ("ch4", "STO-3G"),
             ("hcl", "STO-3G"), ("h2s", "STO-3G")]

CONVENTIONS = {
    "order_spherical": "p: x y z; l >= 2: m = -l..+l (real solid harmonics)",
    "order_cartesian": "lexicographic: xx xy xz yy yz zz, ...",
    "normalisation": "spherical: 1; cartesian l >= 2: radial part only (libcint)",
    "nuclear": "-sum_A Z_A <mu|1/r_A|nu>, point nuclei",
    "dipole": "<mu|r_k|nu>, origin (0,0,0)",
    "eri": "(mu nu|lambda sigma), chemists' notation",
    "packing_symmetric": "lower triangle, row by row: (0,0) (1,0) (1,1) (2,0) ...",
    "packing_eri": "8-fold: ij = i(i+1)/2 + j for i >= j, value index ij(ij+1)/2 + kl for ij >= kl",
}


def tag(spherical: bool) -> str:
    return "sph" if spherical else "cart"


def file_stem(basis_name: str) -> str:
    return basis_name.lower().replace("*", "_st")


def pack_symmetric(matrix: numpy.ndarray) -> list:
    assert numpy.abs(matrix - matrix.T).max() < 1e-12
    return matrix[numpy.tril_indices(matrix.shape[0])].tolist()


def one_electron_matrices(mol) -> dict:
    with mol.with_common_orig((0.0, 0.0, 0.0)):
        dipole = mol.intor("int1e_r", comp=3)
    return {
        "overlap": pack_symmetric(mol.intor("int1e_ovlp")),
        "kinetic": pack_symmetric(mol.intor("int1e_kin")),
        "nuclear": pack_symmetric(mol.intor("int1e_nuc")),
        "dipole_x": pack_symmetric(dipole[0]),
        "dipole_y": pack_symmetric(dipole[1]),
        "dipole_z": pack_symmetric(dipole[2]),
    }


def molecule_record(symbols, coords, charge=0) -> dict:
    return {"charge": charge, "symbols": list(symbols), "geometry_bohr": [list(c) for c in coords]}


def nuclear_repulsion(symbols, coords) -> float:
    """Plain double loop, independent of PySCF."""
    charges = [common.ELEMENTS.index(s) + 1 for s in symbols]
    r = numpy.array(coords, dtype=float)
    return float(sum(charges[a] * charges[b] / numpy.linalg.norm(r[a] - r[b])
                     for a in range(len(charges)) for b in range(a)))


# --------------------------------------------------------------------------

def gen_ints1e() -> None:
    cases = []
    for stem in common.SMALL_MOLECULES:
        for basis_name in HF_BASIS_SETS:
            cases.append((stem, basis_name, common.basis_is_spherical(basis_name)))
    cases += [c for c in EXTRA_1E if c not in cases]

    for stem, basis_name, spherical in cases:
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        if mol.nao > MAX_FUNCTIONS_1E:
            continue
        common.write_json(OUT / f"ints1e_{stem}_{file_stem(basis_name)}_{tag(spherical)}.json", {
            "meta": common.meta("gen_m1.py", conventions=CONVENTIONS),
            "molecule": molecule_record(symbols, coords),
            "basis": {"name": basis_name, "spherical": spherical},
            "n_functions": int(mol.nao),
            "function_labels": [label.strip() for label in mol.ao_labels()],
            "matrices": one_electron_matrices(mol),
        })


def gen_eri() -> None:
    for stem, basis_name in ERI_CASES:
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        n = int(mol.nao)
        packed = mol.intor("int2e", aosym="s8")
        # Make sure the documented packing is the one PySCF uses.
        full = mol.intor("int2e", aosym="s1").reshape(n, n, n, n)
        pair = lambda i, j: i * (i + 1) // 2 + j if i >= j else j * (j + 1) // 2 + i
        rng = numpy.random.default_rng(0)
        for i, j, k, l in rng.integers(0, n, size=(200, 4)):
            ij, kl = pair(i, j), pair(k, l)
            index = ij * (ij + 1) // 2 + kl if ij >= kl else kl * (kl + 1) // 2 + ij
            assert packed[index] == full[i, j, k, l] or abs(packed[index] - full[i, j, k, l]) < 1e-14
        n_pair = n * (n + 1) // 2
        assert packed.size == n_pair * (n_pair + 1) // 2
        common.write_json(OUT / f"eri_{stem}_{file_stem(basis_name)}_{tag(spherical)}.json", {
            "meta": common.meta("gen_m1.py", conventions=CONVENTIONS),
            "molecule": molecule_record(symbols, coords),
            "basis": {"name": basis_name, "spherical": spherical},
            "n_functions": n,
            "eri_packed": packed.tolist(),
        })


# --------------------------------------------------------------------------
# Artificial basis with angular momentum up to l = 6 (i functions).
#
# Three atoms at general positions (no symmetry element, so no integral
# vanishes by symmetry). Exponents and coefficients are arbitrary test values.
# Every pair (l_a, l_b) with l up to 6 occurs between two different centres.

HIGHL_SYMBOLS = ["O", "H", "H"]
HIGHL_COORDS = [[0.10, -0.30, 0.20], [1.30, 0.90, -0.45], [-0.85, 1.15, 1.25]]  # Bohr

HIGHL_BASIS_SPH = """
BASIS "ao basis" SPHERICAL
O    S
      1.30   0.6
      0.45   0.5
O    P
      1.10   1.0
O    D
      1.60   0.4
      0.70   0.7
O    F
      0.90   1.0
O    G
      1.20   0.3
      0.60   0.8
O    H
      0.75   1.0
O    I
      1.00   0.5
      0.55   0.6
H    S
      0.80   1.0
H    F
      0.70   1.0
H    G
      0.65   1.0
H    H
      0.85   0.6
      0.50   0.5
H    I
      0.60   1.0
END
"""

# For the cartesian case one centre carries all angular momenta; the other
# two carry s and p only (identical in both function types).
HIGHL_BASIS_CART = """
BASIS "ao basis" CARTESIAN
O    S
      1.30   0.6
      0.45   0.5
O    P
      1.10   1.0
O    D
      1.60   0.4
      0.70   0.7
O    F
      0.90   1.0
O    G
      1.20   0.3
      0.60   0.8
O    H
      0.75   1.0
O    I
      1.00   0.5
      0.55   0.6
H    S
      0.80   1.0
H    P
      0.70   1.0
END
"""


def element_block(text: str, symbol: str) -> str:
    """Lines of one element from NWChem basis text: headers `<symbol> <type>`
    and the data lines that follow them."""
    selected, active = [], False
    for line in text.splitlines():
        fields = line.split()
        if not fields or fields[0].upper() in ("BASIS", "END"):
            continue
        if fields[0][0].isalpha():
            active = fields[0].upper() == symbol.upper()
        if active:
            selected.append(line)
    assert selected, symbol
    return "\n".join(selected)


def highl_mol(text: str, cart: bool, second_h_minimal: bool):
    basis = {"O": parse_nwchem.parse(element_block(text, "O"), optimize=False),
             "H": parse_nwchem.parse(element_block(text, "H"), optimize=False)}
    if second_h_minimal:
        # Third atom: s shell only, to keep the matrices small.
        basis = {"O": basis["O"], "H1": basis["H"], "H2": [basis["H"][0]]}
        symbols = ["O", "H1", "H2"]
    else:
        symbols = HIGHL_SYMBOLS
    return common.build_mol(symbols, HIGHL_COORDS, basis, cart=cart)


def shell_index(mol, atom: int, l: int) -> int:
    """Index of the (single) shell with angular momentum l on an atom."""
    ids = [i for i in mol.atom_shell_ids(atom) if mol.bas_angular(i) == l]
    assert len(ids) == 1
    return ids[0]


def quartet(mol, shells) -> dict:
    i, j, k, l = shells
    block = mol.intor("int2e", shls_slice=(i, i + 1, j, j + 1, k, k + 1, l, l + 1))
    return {"shells": [int(s) for s in shells], "shape": list(block.shape),
            "values": block.reshape(-1).tolist()}  # C order: last index runs fastest


def gen_highl() -> None:
    for spherical, text in ((True, HIGHL_BASIS_SPH), (False, HIGHL_BASIS_CART)):
        mol = highl_mol(text, cart=not spherical, second_h_minimal=spherical)
        O = lambda l: shell_index(mol, 0, l)
        H1 = lambda l: shell_index(mol, 1, l)
        H2 = lambda l: shell_index(mol, 2, l)
        if spherical:
            quartets = [
                (O(6), H1(6), H2(0), H1(0)),   # (i i|s s), three centres
                (O(6), H1(6), O(6), H2(0)),    # (i i|i s)
                (O(5), H1(5), O(4), H1(3)),    # (h h|g f)
                (O(3), H1(4), H2(0), H1(5)),   # (f g|s h)
                (O(4), O(3), O(2), O(1)),      # one centre
                (O(6), O(6), O(0), O(0)),      # one centre, (i i|s s)
                (O(6), O(5), H1(0), H2(0)),    # (i h|s s)
                (O(2), H1(6), O(1), H1(4)),    # (d i|p g)
            ]
        else:
            quartets = [
                (O(6), H1(0), H2(1), O(2)),    # (i s|p d)
                (O(4), O(3), O(2), H2(1)),     # (g f|d p)
                (O(6), O(0), O(5), O(0)),      # one centre, (i s|h s)
                (O(6), O(6), H1(0), H1(0)),    # (i i|s s)
                (O(5), H1(1), O(3), H2(1)),    # (h p|f p)
            ]
        # The basis as elliDFT has to see it: per atom, in NWChem format.
        atom_basis = []
        for ia in range(mol.natm):
            atom_basis.append([{k: s[k] for k in ("l", "exponents", "coefficients")}
                               for s in common.expanded_shells(mol, atom=ia)])
        common.write_json(OUT / f"highl_{tag(spherical)}.json", {
            "meta": common.meta("gen_m1.py", conventions=CONVENTIONS,
                                note="artificial basis, arbitrary exponents, l = 0..6"),
            "molecule": molecule_record(HIGHL_SYMBOLS, HIGHL_COORDS),
            "basis": {"name": "highl", "spherical": spherical, "shells_per_atom": atom_basis},
            "n_functions": int(mol.nao),
            "n_shells": int(mol.nbas),
            "matrices": one_electron_matrices(mol),
            "eri_shell_quartets": [quartet(mol, q) for q in quartets],
        })


# --------------------------------------------------------------------------

def gen_hf() -> None:
    settings = {
        "method": "RHF, closed shell, exact four-centre integrals (no density fitting)",
        "conv_tol": 1e-12,
        "conv_tol_grad": 1e-9,
        "direct_scf_tol": 1e-14,
        "max_cycle": 200,
        "init_guess": "minao",
        "diis": "PySCF default (CDIIS)",
        "symmetry": False,
        "linear_dependency_removal": "none (all functions kept; see n_orbitals)",
    }
    records = []
    for stem in common.SMALL_MOLECULES:
        symbols, coords, _ = common.read_xyz_bohr(stem)
        for basis_name in HF_BASIS_SETS:
            spherical = common.basis_is_spherical(basis_name)
            mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
            mf = scf.RHF(mol)
            mf.conv_tol = settings["conv_tol"]
            mf.conv_tol_grad = settings["conv_tol_grad"]
            mf.direct_scf_tol = settings["direct_scf_tol"]
            mf.max_cycle = settings["max_cycle"]
            mf.init_guess = settings["init_guess"]
            mf.kernel()
            assert mf.converged, (stem, basis_name)
            assert not hasattr(mf, "with_df"), "density fitting must be off"

            dm = mf.make_rdm1()
            hcore = mf.get_hcore()
            overlap = mf.get_ovlp()
            e_one = float(numpy.einsum("ij,ji->", hcore, dm))
            e_nuc = nuclear_repulsion(symbols, coords)
            assert abs(e_nuc - mf.energy_nuc()) < 1e-11
            n_occ = mol.nelectron // 2
            record = {
                "molecule": stem,
                "basis": basis_name,
                "spherical": spherical,
                "charge": 0,
                "symbols": symbols,
                "geometry_bohr": coords,
                "n_electrons": int(mol.nelectron),
                "n_functions": int(mol.nao),
                "n_orbitals": int(mf.mo_energy.size),
                "overlap_min_eigenvalue": float(numpy.linalg.eigvalsh(overlap)[0]),
                "energy_total": float(mf.e_tot),
                "energy_nuclear_repulsion": e_nuc,
                "energy_one_electron": e_one,
                "energy_two_electron": float(mf.e_tot - e_nuc - e_one),
                "dipole_au": [float(x) for x in mf.dip_moment(unit="AU", verbose=0)],
                "orbital_energies": [float(e) for e in mf.mo_energy],
                "homo": float(mf.mo_energy[n_occ - 1]),
                "lumo": float(mf.mo_energy[n_occ]) if n_occ < mf.mo_energy.size else None,
                "trace_DS": float(numpy.einsum("ij,ji->", dm, overlap)),
                "scf_cycles": int(mf.cycles),
            }
            records.append(record)
            print(f"  {stem:8s} {basis_name:12s} nao {mol.nao:4d}  E = {mf.e_tot:.12f}  "
                  f"cycles {mf.cycles}")
    common.write_json(OUT / "hf.json", {
        "meta": common.meta("gen_m1.py", **settings),
        "results": records,
    })


def main() -> None:
    parts = {"ints1e": gen_ints1e, "eri": gen_eri, "highl": gen_highl, "hf": gen_hf}
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", choices=sorted(parts), help="generate one part only")
    args = parser.parse_args()
    for name, function in parts.items():
        if args.only in (None, name):
            function()


if __name__ == "__main__":
    main()
