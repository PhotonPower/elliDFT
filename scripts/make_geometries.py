#!/usr/bin/env python3
"""Generate the test geometries of elliDFT (data/geometries/*.xyz, in Angstrom).

Two groups:

1. Small molecules, built here from nominal internal coordinates (bond
   lengths in Angstrom, angles in degrees). The values are round,
   near-experimental numbers chosen as *test inputs*. They are not literature
   values and must not be cited as such. elliDFT and the reference code get
   exactly the same coordinates, so their accuracy as molecular structures
   does not enter any acceptance criterion.
   Exception: H2 is placed at exactly 1.4 Bohr (anchor of the specification).

2. Repeat units and oligomers of polystyrene and PMMA (hydrogen terminated),
   generated with RDKit: extended chains with an all-trans backbone along z,
   syndiotactic, MMFF94 with the backbone torsions held at 180 degrees. All
   oligomers of a series share one orientation, so that the polarisability
   increments refer to one chain axis. Isotactic chains with 1 to 3 units
   (suffix _iso) serve as a control. This follows the decision of the client
   of 2026-11-14 (docs/Entscheidungen.md). MMFF94 is a force field; the
   geometries are not optimised at the level of the calculations.

Usage:  python3 scripts/make_geometries.py            # small molecules only
        python3 scripts/make_geometries.py --oligomers  # also group 2 (needs RDKit)
"""

import argparse
import math
import pathlib
import sys

BOHR_IN_ANGSTROM = 0.529177210903  # CODATA 2018, as fixed in the specification

OUT_DIR = pathlib.Path(__file__).resolve().parent.parent / "data" / "geometries"


def write_xyz(name: str, comment: str, atoms: list[tuple[str, float, float, float]]) -> None:
    lines = [str(len(atoms)), comment]
    for symbol, x, y, z in atoms:
        # +0.0 turns -0.0 into 0.0
        lines.append(f"{symbol:<2s} {x + 0.0:21.15f} {y + 0.0:21.15f} {z + 0.0:21.15f}")
    (OUT_DIR / f"{name}.xyz").write_text("\n".join(lines) + "\n")
    print(f"wrote {name}.xyz ({len(atoms)} atoms)")


# --------------------------------------------------------------------------
# Group 1: small molecules from nominal internal coordinates
# --------------------------------------------------------------------------

def small_molecules() -> dict[str, tuple[str, list]]:
    mols = {}
    nominal = "nominal test geometry (scripts/make_geometries.py), not a literature value"

    r = 1.4 * BOHR_IN_ANGSTROM
    mols["h2"] = ("H2, bond length exactly 1.4 Bohr (anchor of the specification)",
                  [("H", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, r)])

    for atom in ("He", "Ne", "Ar"):
        mols[atom.lower()] = (f"{atom} atom at the origin", [(atom, 0.0, 0.0, 0.0)])

    def bent(center: str, r: float, angle_deg: float):
        """AH2 in the yz plane, C2 axis along z, heavy atom at the origin."""
        half = math.radians(angle_deg) / 2.0
        return [(center, 0.0, 0.0, 0.0),
                ("H", 0.0, r * math.sin(half), r * math.cos(half)),
                ("H", 0.0, -r * math.sin(half), r * math.cos(half))]

    mols["h2o"] = (f"H2O, r(OH) = 0.9572, HOH = 104.52; {nominal}", bent("O", 0.9572, 104.52))
    mols["h2s"] = (f"H2S, r(SH) = 1.3356, HSH = 92.12; {nominal}", bent("S", 1.3356, 92.12))

    # NH3: C3 axis along z, N at the origin. For three equivalent bonds with
    # mutual angle t the polar angle b obeys sin(b) = 2 sin(t/2) / sqrt(3).
    r, t = 1.012, math.radians(106.67)
    sin_b = 2.0 * math.sin(t / 2.0) / math.sqrt(3.0)
    cos_b = -math.sqrt(1.0 - sin_b * sin_b)  # hydrogens below the nitrogen
    nh3 = [("N", 0.0, 0.0, 0.0)]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        nh3.append(("H", r * sin_b * math.cos(phi), r * sin_b * math.sin(phi), r * cos_b))
    mols["nh3"] = (f"NH3, r(NH) = 1.012, HNH = 106.67; {nominal}", nh3)

    # CH4: tetrahedral, C at the origin, hydrogens on alternating cube corners.
    d = 1.087 / math.sqrt(3.0)
    mols["ch4"] = (f"CH4, r(CH) = 1.087, tetrahedral; {nominal}",
                   [("C", 0.0, 0.0, 0.0), ("H", d, d, d), ("H", d, -d, -d),
                    ("H", -d, d, -d), ("H", -d, -d, d)])

    mols["hf"] = (f"HF, r = 0.9168; {nominal}", [("F", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 0.9168)])
    mols["hcl"] = (f"HCl, r = 1.2746; {nominal}", [("Cl", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 1.2746)])

    # C2H4: planar (yz plane), C=C along z, centred at the origin.
    rcc, rch, a = 1.339, 1.087, math.radians(121.3)
    c2h4 = [("C", 0.0, 0.0, rcc / 2.0), ("C", 0.0, 0.0, -rcc / 2.0)]
    for sz in (1.0, -1.0):
        for sy in (1.0, -1.0):
            c2h4.append(("H", 0.0, sy * rch * math.sin(a), sz * (rcc / 2.0 - rch * math.cos(a))))
    mols["c2h4"] = (f"C2H4, r(CC) = 1.339, r(CH) = 1.087, HCC = 121.3; {nominal}", c2h4)

    # Benzene: regular hexagon in the xy plane, centred at the origin.
    rcc, rch = 1.397, 1.084
    benzene = []
    for symbol, radius in (("C", rcc), ("H", rcc + rch)):
        for k in range(6):
            phi = math.pi * k / 3.0
            benzene.append((symbol, radius * math.cos(phi), radius * math.sin(phi), 0.0))
    mols["benzene"] = (f"C6H6, r(CC) = 1.397, r(CH) = 1.084, D6h; {nominal}", benzene)
    return mols


# --------------------------------------------------------------------------
# Group 2: repeat units and oligomers via RDKit
#
# Decision of the client (2026-11-14): extended chains (all-trans backbone)
# along z, syndiotactic, built with RDKit and MMFF94; isotactic chains with 1
# to 3 units as a control. All oligomers of a series share one orientation:
# a template one unit longer than the longest oligomer is built first; all
# oligomers are sections of it (heavy atoms held at its positions, hydrogens
# relaxed), so that unit k of every oligomer lies
# where unit k of the longest one lies.
# --------------------------------------------------------------------------

EMBED_SEED = 20261114

# Per polymer: atoms of one repeat unit in SMILES order. The first atom is the
# substituted backbone carbon, `side` is the offset of the atom that defines
# the side of the substituent (phenyl ipso carbon, carbonyl carbon).
POLYMERS = {
    "styrene": {
        "unit": "{c}(c1ccccc1)", "stereo": ("[C@H]", "[C@@H]"), "last": "C(c1ccccc1)",
        "atoms_per_unit": 7, "side": 1,
        "monomer_name": "ethylbenzene", "description": "styrene",
    },
    "mma": {
        "unit": "{c}(C)(C(=O)OC)", "stereo": ("[C@]", "[C@@]"), "last": "C(C)(C(=O)OC)",
        "atoms_per_unit": 6, "side": 2,
        "monomer_name": "methyl_isobutyrate", "description": "MMA",
    },
}


def oligomer_smiles(polymer: str, n: int, tacticity: str | None) -> str:
    """H-[CH2-CXY]n-H written as CH3-CXY-CH2-CXY-...-CH2-CHXY.

    With `tacticity` the stereocentres of all units but the last carry
    chirality marks: the same mark for every unit, or alternating marks. Which
    of the two is isotactic is not assumed here but measured on the finished
    geometry (see tacticity_signs). Without tacticity no marks are written.
    """
    spec = POLYMERS[polymer]
    parts = ["C"]
    for k in range(n - 1):
        if tacticity is None:
            centre = "C"
        else:
            centre = spec["stereo"][0] if (tacticity == "same" or k % 2 == 0) else spec["stereo"][1]
        parts.append(spec["unit"].format(c=centre) + "C")
    parts.append(spec["last"])
    return "".join(parts)


def backbone_indices(polymer: str, n: int) -> list[int]:
    """Indices of the backbone carbons CH3, CXY, CH2, CXY, ..., CHXY (heavy atoms come first)."""
    per_unit = POLYMERS[polymer]["atoms_per_unit"]
    indices = [0]
    for k in range(n):
        centre = 1 + k * (per_unit + 1)
        indices.append(centre)
        if k < n - 1:
            indices.append(centre + per_unit)
    return indices


def side_indices(polymer: str, n: int) -> list[tuple[int, int]]:
    """(substituted backbone carbon, atom that defines the side of the substituent) per unit."""
    spec = POLYMERS[polymer]
    return [(1 + k * (spec["atoms_per_unit"] + 1), 1 + k * (spec["atoms_per_unit"] + 1) + spec["side"])
            for k in range(n)]


def constrained_optimisation(mol, backbone: list[int], pre_rotate: bool, hold_heavy_atoms: bool = False) -> None:
    """MMFF94 with all backbone torsions held at 180 degrees.

    With `pre_rotate` the backbone is first turned rigidly into the all-trans
    conformation, so that the force field only has to relax the side groups.
    With `hold_heavy_atoms` every heavy atom is tied to its start position
    (free within 0.05 Angstrom), so that essentially only hydrogens move."""
    from rdkit.Chem import AllChem, rdMolTransforms
    conformer = mol.GetConformer()
    torsions = list(zip(backbone, backbone[1:], backbone[2:], backbone[3:]))
    if pre_rotate:
        for a, b, c, d in torsions:
            rdMolTransforms.SetDihedralDeg(conformer, a, b, c, d, 180.0)
    properties = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    field = AllChem.MMFFGetMoleculeForceField(mol, properties)
    for a, b, c, d in torsions:
        field.MMFFAddTorsionConstraint(a, b, c, d, False, 180.0, 180.0, 1.0e3)
    if hold_heavy_atoms:
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() > 1:
                field.MMFFAddPositionConstraint(atom.GetIdx(), 0.05, 1.0e3)
    for _ in range(50):
        if field.Minimize(maxIts=5000) == 0:
            return
    raise RuntimeError("constrained MMFF94 optimisation not converged")


def chain_orientation(positions, backbone: list[int], first_side: tuple[int, int]):
    """Rotation and shift that put the backbone along z, the zigzag plane into
    xz and the first substituent on the +y side; returns transformed positions."""
    import numpy
    chain = positions[backbone]
    centre = chain.mean(axis=0)
    # Principal axes of the backbone carbons: largest extent = chain axis,
    # smallest = normal of the zigzag plane.
    _, vectors = numpy.linalg.eigh((chain - centre).T @ (chain - centre))
    z, y = vectors[:, 2], vectors[:, 0]
    if (chain[-1] - chain[0]) @ z < 0.0:
        z = -z
    if (positions[first_side[1]] - positions[first_side[0]]) @ y < 0.0:
        y = -y
    x = numpy.cross(y, z)
    return (positions - centre) @ numpy.column_stack([x, y, z])


def superimpose(moving, reference_points, indices):
    """Kabsch: rotates and shifts `moving` so that moving[indices] matches reference_points."""
    import numpy
    a, b = moving[indices], reference_points
    ca, cb = a.mean(axis=0), b.mean(axis=0)
    u, _, vt = numpy.linalg.svd((a - ca).T @ (b - cb))
    d = numpy.sign(numpy.linalg.det(u @ vt))
    rotation = u @ numpy.diag([1.0, 1.0, d]) @ vt
    return (moving - ca) @ rotation + cb


def tacticity_signs(positions, backbone: list[int], sides: list[tuple[int, int]]) -> list[int]:
    """Side of the zigzag plane (+1 or -1) on which the substituent of each unit lies."""
    import numpy
    chain = positions[backbone]
    centre = chain.mean(axis=0)
    _, vectors = numpy.linalg.eigh((chain - centre).T @ (chain - centre))
    normal = vectors[:, 0]
    return [int(numpy.sign((positions[s] - positions[c]) @ normal)) for c, s in sides]


def backbone_torsions(mol, backbone: list[int]) -> list[float]:
    from rdkit.Chem import rdMolTransforms
    conformer = mol.GetConformer()
    return [abs(rdMolTransforms.GetDihedralDeg(conformer, a, b, c, d))
            for a, b, c, d in zip(backbone, backbone[1:], backbone[2:], backbone[3:])]


def build_series(polymer: str, sizes: list[int], pattern: str, label: str, suffix: str) -> None:
    """Builds the oligomers `sizes` of one polymer and tacticity pattern."""
    import numpy
    import rdkit
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    spec = POLYMERS[polymer]
    # A template one unit longer than the longest oligomer wanted: its last
    # unit is no stereocentre (the embedding puts its substituent on either
    # side), so it is never written; every oligomer is a section of the
    # template made of units whose side is fixed by the tacticity.
    n_max = max(sizes) + 1

    # Longest oligomer: stereocentres from the SMILES, then all-trans backbone.
    # The tacticity is read off the finished geometry: in an all-trans chain
    # the substituents of an isotactic polymer lie on the same side of the
    # zigzag plane, those of a syndiotactic polymer alternate. The last unit
    # is not a stereocentre of its own and is left out. The first unit of an
    # MMA oligomer carries two methyl groups and is no stereocentre either;
    # there the embedding decides, so seeds are tried until the sides match.
    backbone = backbone_indices(polymer, n_max)
    longest, reference, seed_used = None, None, None
    for seed in range(EMBED_SEED, EMBED_SEED + 40):
        candidate = Chem.AddHs(Chem.MolFromSmiles(oligomer_smiles(polymer, n_max, pattern)))
        if AllChem.EmbedMolecule(candidate, randomSeed=seed) != 0:
            continue
        constrained_optimisation(candidate, backbone, pre_rotate=True)
        positions = chain_orientation(candidate.GetConformer().GetPositions(), backbone,
                                      side_indices(polymer, n_max)[0])
        signs = tacticity_signs(positions, backbone, side_indices(polymer, n_max))[:-1]
        same = all(s == signs[0] for s in signs)
        alternating = all(a == -b for a, b in zip(signs, signs[1:]))
        if len(signs) < 2 or {"isotactic": same, "syndiotactic": alternating}[label]:
            longest, reference, seed_used = candidate, positions, seed
            break
    if longest is None:
        raise RuntimeError(f"{polymer} {label}: substituent sides do not match for any seed")

    for n in sorted(sizes, reverse=True):
        # A section of the template: the heavy atoms of the n-mer are the
        # first heavy atoms of the template and keep their positions (within
        # 0.05 Angstrom); hydrogens are added and relaxed. The increments
        # alpha(m) - alpha(m-1) then describe the addition of one unit, not a
        # change of conformation.
        bare = Chem.MolFromSmiles(oligomer_smiles(polymer, n, None))
        n_heavy = bare.GetNumAtoms()
        conformer = Chem.Conformer(n_heavy)
        for i in range(n_heavy):
            conformer.SetAtomPosition(i, Point3D(*reference[i]))
        bare.AddConformer(conformer, assignId=True)
        mol = Chem.AddHs(bare, addCoords=True)
        constrained_optimisation(mol, backbone_indices(polymer, n), pre_rotate=False, hold_heavy_atoms=True)
        heavy = list(range(n_heavy))
        positions = superimpose(mol.GetConformer().GetPositions(), reference[heavy], heavy)
        deviation = numpy.abs(positions[heavy] - reference[heavy]).max()
        if deviation > 0.15:
            raise RuntimeError(f"{polymer} {n}: heavy atoms deviate by {deviation:.2f} A from the template")
        # Measured on the written geometry: every unit on the side the tacticity demands.
        written = tacticity_signs(positions, backbone_indices(polymer, n), side_indices(polymer, n)) if n > 1 else [1]
        ok = all(s == written[0] for s in written) if label == "isotactic" else \
            all(a == -b for a, b in zip(written, written[1:]))
        if not ok:
            raise RuntimeError(f"{polymer} {n} {label}: substituent sides {written}")
        torsions = backbone_torsions(mol, backbone_indices(polymer, n))
        if torsions and min(torsions) < 170.0:
            raise RuntimeError(f"{polymer} {n}: backbone torsion {min(torsions):.1f} is not trans")
        atoms = [(a.GetSymbol(), *positions[a.GetIdx()]) for a in mol.GetAtoms()]
        name = (spec["monomer_name"] if n == 1 else f"{polymer}_{n}mer") + suffix
        comment = (f"{spec['description']} oligomer, {n} unit{'s' if n > 1 else ''}, H-terminated, {label}, "
                   f"all-trans backbone along z, zigzag plane xz; RDKit {rdkit.__version__}, MMFF94 with backbone "
                   f"torsions held at 180 deg; section of a {n_max}-mer template (embedding seed {seed_used}), "
                   f"heavy atoms held (largest deviation {deviation:.2f} A), hydrogens relaxed")
        write_xyz(name, comment, atoms)


def oligomers() -> None:
    # Which chirality pattern is syndiotactic is determined on the geometry:
    # try the alternating pattern first and fall back to the other one.
    for polymer, sizes in (("styrene", list(range(1, 9))), ("mma", list(range(1, 5)))):
        built = {}
        for label, suffix, wanted_sizes in (("syndiotactic", "", sizes), ("isotactic", "_iso", [1, 2, 3])):
            for pattern in ("alternating", "same"):
                if pattern in built.values():
                    continue
                try:
                    build_series(polymer, wanted_sizes, pattern, label, suffix)
                    built[label] = pattern
                    break
                except RuntimeError as error:
                    if "do not match" not in str(error):
                        raise
            if label not in built:
                raise RuntimeError(f"{polymer}: no chirality pattern gives a {label} chain")
        print(f"{polymer}: syndiotactic = {built['syndiotactic']} marks, isotactic = {built['isotactic']} marks")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--oligomers", action="store_true",
                        help="also generate the oligomer geometries (needs RDKit)")
    args = parser.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, (comment, atoms) in small_molecules().items():
        write_xyz(name, comment, atoms)
    if args.oligomers:
        oligomers()
    return 0


if __name__ == "__main__":
    sys.exit(main())
