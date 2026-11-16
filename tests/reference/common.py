"""Shared helpers of the reference scripts (PySCF side).

Conventions, set explicitly so that both programs compute the same thing:

* Geometries are read from data/geometries/*.xyz (Angstrom), converted with
  the CODATA 2018 value fixed in the specification and passed to PySCF in
  Bohr. PySCF's own conversion constant is never used.
* Basis sets are read from the files in data/basis (the same files elliDFT
  reads), never from PySCF's built-in library. Parsing uses optimize=False, so
  PySCF keeps the file order within one angular momentum.
* Spherical or cartesian functions are always passed explicitly (`cart`).
* No point group symmetry, no density fitting unless a script says so.
"""

import json
import pathlib
import platform
import re

import numpy
import pyscf
import scipy
from pyscf import gto
from pyscf.gto.basis import parse_nwchem

BOHR_IN_ANGSTROM = 0.529177210903  # CODATA 2018, specification section 5

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
BASIS_DIR = ROOT / "data" / "basis"
GEOMETRY_DIR = ROOT / "data" / "geometries"
DATA_DIR = pathlib.Path(__file__).resolve().parent / "data"

ORBITAL_BASIS_SETS = ["STO-3G", "6-31G*", "def2-SVP", "def2-TZVP", "def2-SVPD", "def2-TZVPD",
                      "cc-pVDZ", "aug-cc-pVDZ", "aug-cc-pVTZ"]
AUXILIARY_BASIS_SETS = ["def2-universal-JKFIT", "def2-universal-JFIT"]
ALL_BASIS_SETS = ORBITAL_BASIS_SETS + AUXILIARY_BASIS_SETS

ELEMENTS = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]

# Test molecules of the specification, section 7 (file stems in data/geometries).
SMALL_MOLECULES = ["h2", "he", "h2o", "nh3", "ch4", "hf", "hcl", "h2s", "ne", "ar",
                   "c2h4", "benzene"]
REPEAT_UNITS = ["ethylbenzene", "methyl_isobutyrate"]
OLIGOMERS = [f"styrene_{n}mer" for n in range(2, 9)] + [f"mma_{n}mer" for n in range(2, 5)]
# Isotactic control series, 1 to 3 units (decision of the client, 2026-11-14).
ISOTACTIC_CONTROLS = ["ethylbenzene_iso", "styrene_2mer_iso", "styrene_3mer_iso",
                      "methyl_isobutyrate_iso", "mma_2mer_iso", "mma_3mer_iso"]
ALL_MOLECULES = SMALL_MOLECULES + REPEAT_UNITS + OLIGOMERS + ISOTACTIC_CONTROLS


def basis_file(basis_name: str) -> pathlib.Path:
    """Same file name convention as src/basis and scripts/fetch_basis.py."""
    return BASIS_DIR / (basis_name.lower().replace("*", "_st") + ".nw")


def basis_is_spherical(basis_name: str) -> bool:
    """Function type declared in the header of the basis set file."""
    for line in basis_file(basis_name).read_text().splitlines():
        if line.strip().upper().startswith("BASIS"):
            return "SPHERICAL" in line.upper()
    raise RuntimeError(f"{basis_name}: no BASIS header")


def read_xyz_bohr(stem: str):
    """Returns (symbols, coordinates in Bohr, comment line)."""
    lines = (GEOMETRY_DIR / f"{stem}.xyz").read_text().splitlines()
    n = int(lines[0])
    symbols, coords = [], []
    for line in lines[2:2 + n]:
        fields = line.split()
        symbols.append(fields[0])
        coords.append([float(x) / BOHR_IN_ANGSTROM for x in fields[1:4]])
    return symbols, coords, lines[1].strip()


def pyscf_basis(basis_name: str, symbols) -> dict:
    """Basis in PySCF's internal format, read from our file, file order kept."""
    path = str(basis_file(basis_name))
    return {s: parse_nwchem.load(path, s, optimize=False) for s in sorted(set(symbols))}


def build_mol(symbols, coords_bohr, basis, cart: bool, charge: int = 0) -> gto.Mole:
    """`basis` is a basis set name (read from data/basis) or a PySCF basis dict."""
    mol = gto.Mole()
    mol.atom = [[s, tuple(xyz)] for s, xyz in zip(symbols, coords_bohr)]
    mol.unit = "Bohr"
    mol.basis = pyscf_basis(basis, symbols) if isinstance(basis, str) else basis
    mol.cart = cart
    mol.charge = charge
    mol.spin = None  # lowest spin; all molecular test systems are closed shell
    mol.symmetry = False
    mol.verbose = 0
    mol.build()
    # Guard: PySCF must have taken the coordinates unchanged.
    assert numpy.array_equal(mol.atom_coords(), numpy.array(coords_bohr, dtype=float).reshape(-1, 3))
    return mol


def expanded_shells(mol: gto.Mole, atom: int = None):
    """Segmented view of the basis as PySCF holds it, in AO order.

    One entry per contraction: {"atom", "l", "exponents", "coefficients"}.
    Coefficients are the raw values of the basis set file (mol._basis);
    primitives with zero coefficient in a contraction are left out. This is
    the representation elliDFT uses (src/basis/basis_library.hpp).
    """
    out = []
    for ia in range(mol.natm):
        if atom is not None and ia != atom:
            continue
        raw = mol._basis[mol.atom_symbol(ia)]
        shell_ids = mol.atom_shell_ids(ia)
        assert len(raw) == len(shell_ids)
        for shell_id, entry in zip(shell_ids, raw):
            l = entry[0]
            assert not isinstance(entry[1], int), "kappa entries are not expected"
            table = numpy.array(entry[1:], dtype=float)
            # Consistency with what the integral library sees.
            assert l == mol.bas_angular(shell_id)
            assert table.shape[1] - 1 == mol.bas_nctr(shell_id)
            assert numpy.array_equal(table[:, 0], mol.bas_exp(shell_id))
            for col in range(1, table.shape[1]):
                keep = table[:, col] != 0.0
                out.append({"atom": ia, "l": int(l),
                            "exponents": table[keep, 0].tolist(),
                            "coefficients": table[keep, col].tolist()})
    return out


def meta(script: str, **settings) -> dict:
    return {
        "generated_by": f"tests/reference/{script}",
        "program": "PySCF",
        "pyscf_version": pyscf.__version__,
        "numpy_version": numpy.__version__,
        "scipy_version": scipy.__version__,
        "python_version": platform.python_version(),
        "bohr_in_angstrom": BOHR_IN_ANGSTROM,
        "settings": settings,
    }


def _compact_arrays(text: str) -> str:
    """Puts innermost arrays of numbers on one line (smaller, readable diffs)."""
    return re.sub(r"\[\s+([^\[\]{}]*?)\s+\]",
                  lambda m: "[" + re.sub(r",\s+", ", ", m.group(1)) + "]", text)


def write_json(path: pathlib.Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_compact_arrays(json.dumps(obj, indent=1)) + "\n")
    print(f"wrote {path.relative_to(ROOT)}  ({path.stat().st_size / 1024:.0f} kB)")
