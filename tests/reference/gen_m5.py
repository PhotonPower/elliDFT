#!/usr/bin/env python3
"""Reference data for milestone M5: linear response (PySCF).

Writes data/m5/response.json with, per molecule and functional:
  - the ten lowest singlet excitations in full linear response (RPA) and in
    the Tamm-Dancoff approximation (TDA): energies and oscillator strengths
    in the length gauge,
  - the static polarisability tensor.

How the numbers are obtained. PySCF builds the matrices A and B of the
response equations explicitly (pyscf.tdscf.rhf.get_ab: orbital energy
differences, Coulomb and exchange integrals including the range-separated
part, exchange-correlation kernel on the grid). They are diagonalised here
with dense linear algebra:

  RPA:  (A-B)^(1/2) (A+B) (A-B)^(1/2) y = omega^2 y,  X+Y = (A-B)^(1/2) y / sqrt(omega)
  TDA:  A X = omega X
  t = sqrt(2) mu^T (X+Y),  f = 2/3 omega |t|^2,  alpha = 4 mu^T (A+B)^-1 mu

This avoids PySCF's iterative solver for the full response problem, which
does not converge reliably (H2O with PBE0: some states; HF with PBE0: none,
with energies off by tenths of a Hartree), and the many SCF runs of a
finite-field polarisability. Cross-checks for one small case (HF molecule,
PBE0), asserted in this script:
  - PySCF's Tamm-Dancoff solver gives the same energies and oscillator
    strengths as the dense diagonalisation,
  - finite fields in PySCF give the same polarisability.
The state of PySCF's full-response solver for that case is recorded in the
output file.

Grid (99, 590), unpruned, as for the energies. Functionals by libxc names.

Usage: python tests/reference/gen_m5.py
"""

import numpy
import scipy.linalg
from pyscf import dft
from pyscf.tdscf import rhf as td_rhf

import common
from gen_m2 import make_grid

OUT = common.DATA_DIR / "m5"
KS_GRID = (99, 590)
FUNCTIONALS = {"PBE": "GGA_X_PBE,GGA_C_PBE", "PBE0": "HYB_GGA_XC_PBEH", "CAM-B3LYP": "HYB_GGA_XC_CAM_B3LYP"}
CASES = [("hf", "def2-SVP"), ("h2o", "def2-SVP"), ("c2h4", "def2-SVP")]
CROSS_CHECK = ("hf", "PBE0")
N_STATES = 10
FIELD = 2e-3


def run_scf(mol, xc, field=None):
    mf = dft.RKS(mol)
    mf.xc = xc
    mf.grids = make_grid(mol, *KS_GRID)
    mf.small_rho_cutoff = 0.0
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.direct_scf_tol = 1e-14
    mf.max_cycle = 200
    if field is not None:
        h0 = mf.get_hcore()
        with mol.with_common_orig((0.0, 0.0, 0.0)):
            r = mol.intor("int1e_r", comp=3)
        h = h0 + numpy.einsum("x,xij->ij", numpy.asarray(field), r)  # electrons: potential +F.r
        mf.get_hcore = lambda *args, **kwargs: h
    mf.kernel()
    assert mf.converged
    return mf


def response_from_ab(mf) -> dict:
    mol = mf.mol
    a, b = td_rhf.get_ab(mf)
    nocc, nvir = a.shape[:2]
    n = nocc * nvir
    A, B = a.reshape(n, n), b.reshape(n, n)
    assert numpy.abs(A - A.T).max() < 1e-8 and numpy.abs(B - B.T).max() < 1e-8
    with mol.with_common_orig((0.0, 0.0, 0.0)):
        r = mol.intor("int1e_r", comp=3)
    C = mf.mo_coeff
    mu = numpy.array([(C[:, :nocc].T @ r[k] @ C[:, nocc:]).reshape(n) for k in range(3)])  # 3 x n

    def strengths(omega, z):
        t = numpy.sqrt(2.0) * mu @ z                         # 3 x states
        return 2.0 / 3.0 * omega * (t * t).sum(axis=0)

    # Full response.
    root = scipy.linalg.sqrtm(A - B).real
    w2, y = numpy.linalg.eigh(root @ (A + B) @ root)
    assert w2.min() > 0.0
    omega = numpy.sqrt(w2)
    z = root @ y / numpy.sqrt(omega)
    rpa = {"energies": omega[:N_STATES].tolist(), "oscillator_strengths": strengths(omega, z)[:N_STATES].tolist()}

    # Tamm-Dancoff.
    w, x = numpy.linalg.eigh(A)
    tda = {"energies": w[:N_STATES].tolist(), "oscillator_strengths": strengths(w, x)[:N_STATES].tolist()}

    alpha = 4.0 * mu @ numpy.linalg.solve(A + B, mu.T)
    return {"rpa": rpa, "tda": tda, "polarizability_static": alpha.tolist(), "n_pairs": int(n)}


def grouped(values, energies):
    """Sums over groups of (nearly) degenerate states."""
    sums, start = [], 0
    for i in range(1, len(energies) + 1):
        if i == len(energies) or energies[i] - energies[i - 1] > 1e-5:
            sums.append(sum(values[start:i]))
            start = i
    return numpy.array(sums)


def cross_check(mf, mol, xc, result) -> dict:
    """PySCF's own TDDFT/TDA solvers and finite fields against the numbers above."""
    report = {}
    # Tamm-Dancoff: PySCF's Davidson solver converges and must reproduce the
    # dense diagonalisation, energies and oscillator strengths.
    td = mf.TDA()
    td.singlet = True
    td.nstates = N_STATES
    td.kernel()
    assert all(td.converged)
    e_ref = numpy.array(result["tda"]["energies"])
    de = numpy.abs(numpy.array(td.e) - e_ref).max()
    df = numpy.abs(grouped(td.oscillator_strength(gauge="length"), e_ref) -
                   grouped(result["tda"]["oscillator_strengths"], e_ref)).max()
    assert de < 1e-8 and df < 1e-6, (de, df)
    report["tda_solver"] = {"max_energy_difference": float(de), "max_strength_difference": float(df)}

    # Full response: PySCF's iterative solver is recorded but not relied on.
    # For this case (degenerate pi states) it does not converge at all.
    td = mf.TDDFT()
    td.singlet = True
    td.nstates = N_STATES
    td.kernel()
    converged = numpy.array(td.converged, dtype=bool)
    report["rpa_solver"] = {
        "states_converged_in_pyscf": int(converged.sum()),
        "max_energy_difference_of_converged_states":
            float(numpy.abs(numpy.array(td.e) - numpy.array(result["rpa"]["energies"]))[converged].max())
            if converged.any() else None,
    }
    alpha = numpy.zeros((3, 3))
    for l in range(3):
        def dipole(strength):
            field = [0.0, 0.0, 0.0]
            field[l] = strength
            return run_scf(mol, xc, field).dip_moment(unit="AU", verbose=0)
        d1 = (dipole(FIELD) - dipole(-FIELD)) / (2 * FIELD)
        d2 = (dipole(2 * FIELD) - dipole(-2 * FIELD)) / (4 * FIELD)
        alpha[:, l] = (4 * d1 - d2) / 3
    exact = numpy.array(result["polarizability_static"])
    deviation = numpy.abs(alpha - exact).max() / numpy.abs(exact).max()
    assert deviation < 1e-5, deviation
    report["finite_field_polarizability_relative_difference"] = float(deviation)
    return report


def main() -> None:
    settings = {
        "functionals": FUNCTIONALS,
        "grid": {"n_radial": KS_GRID[0], "n_angular": KS_GRID[1]},
        "small_rho_cutoff": 0.0,
        "n_states": N_STATES,
        "source": "A and B from pyscf.tdscf.rhf.get_ab, dense diagonalisation; singlets; length gauge",
        "cross_check_case": list(CROSS_CHECK),
    }
    records = []
    for stem, basis_name in CASES:
        spherical = common.basis_is_spherical(basis_name)
        symbols, coords, _ = common.read_xyz_bohr(stem)
        mol = common.build_mol(symbols, coords, basis_name, cart=not spherical)
        for label, xc in FUNCTIONALS.items():
            mf = run_scf(mol, xc)
            result = response_from_ab(mf)
            record = {
                "molecule": stem, "basis": basis_name, "spherical": spherical, "functional": label,
                "charge": 0, "symbols": symbols, "geometry_bohr": coords,
                "n_functions": int(mol.nao), "energy_total": float(mf.e_tot), **result,
            }
            if (stem, label) == CROSS_CHECK:
                settings["cross_check"] = cross_check(mf, mol, xc, result)
            records.append(record)
            print(f"  {stem:8s} {label:10s} lowest RPA {result['rpa']['energies'][0]:.6f}  "
                  f"TDA {result['tda']['energies'][0]:.6f}  alpha_zz {result['polarizability_static'][2][2]:.5f}",
                  flush=True)
    common.write_json(OUT / "response.json", {"meta": common.meta("gen_m5.py", **settings), "results": records})


if __name__ == "__main__":
    main()
