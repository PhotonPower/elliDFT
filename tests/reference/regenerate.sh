#!/usr/bin/env bash
# Regenerates the reference data under tests/reference/data with one command.
#
#   tests/reference/regenerate.sh m0     basis structure (PySCF)
#   tests/reference/regenerate.sh m1     Boys function, integrals, Hartree-Fock (PySCF)
#   tests/reference/regenerate.sh m2     grids, functionals, Kohn-Sham energies (PySCF)
#   tests/reference/regenerate.sh m3     density fitting: integrals and energies (PySCF)
#   tests/reference/regenerate.sh m4     hybrid functionals, erf-attenuated integrals (PySCF)
#   tests/reference/regenerate.sh m5     linear response: excitations, polarisability (PySCF)
#   tests/reference/regenerate.sh m6     measured refractive indices of polystyrene and PMMA (no PySCF)
#   tests/reference/regenerate.sh all
#
# Python with the packages of tests/reference/requirements.txt is taken from
# $PYTHON, or from .venv in the repository root. To set it up:
#   python3 -m venv .venv && .venv/bin/pip install -r tests/reference/requirements.txt
#
# The tests themselves do not need Python or PySCF; they read the stored JSON files.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "$here/../.." && pwd)"
python="${PYTHON:-$root/.venv/bin/python}"

if ! "$python" -c "import pyscf" >/dev/null 2>&1; then
  echo "error: '$python' cannot import pyscf; see the header of this script" >&2
  exit 1
fi

what="${1:-}"
case "$what" in
  m0)
    "$python" "$here/gen_m0.py"
    ;;
  m1)
    "$python" "$here/gen_boys.py"
    "$python" "$here/gen_m1.py"
    ;;
  m2)
    "$python" "$here/gen_m2.py"
    ;;
  m3)
    "$python" "$here/gen_m3.py"
    ;;
  m4)
    "$python" "$here/gen_m4.py"
    ;;
  m5)
    "$python" "$here/gen_m5.py"
    ;;
  m6)
    "$python" "$here/gen_m6.py"
    ;;
  all)
    "$0" m0
    "$0" m1
    "$0" m2
    "$0" m3
    "$0" m4
    "$0" m5
    "$0" m6
    ;;
  *)
    echo "usage: $0 m0|m1|m2|m3|m4|m5|m6|all" >&2
    exit 2
    ;;
esac
