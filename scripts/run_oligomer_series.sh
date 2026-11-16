#!/usr/bin/env bash
# Runs the response calculation for an oligomer series and derives n and k.
#
#   scripts/run_oligomer_series.sh <name> <density g/cm^3> <molar mass g/mol> <output dir> \
#       <1 unit>.xyz <2 units>.xyz [<3 units>.xyz ...]
#
# Settings through environment variables (defaults follow the specification, M6):
#   FUNCTIONAL=CAM-B3LYP  BASIS=def2-SVPD  AUX_BASIS=def2-universal-JKFIT  GRID=75,302
#   WAVELENGTHS=400,450,500,550,589,600,650,700,750,800,850,900,950,1000   (nm)
#   GAMMA=0  STATES=0  ORDER_PARAMETER=0  BUILD_DIR=build
#
# Example (polystyrene, guide density of the specification):
#   scripts/run_oligomer_series.sh polystyrene 1.05 104.15 results/ps \
#       data/geometries/ethylbenzene.xyz data/geometries/styrene_2mer.xyz \
#       data/geometries/styrene_3mer.xyz data/geometries/styrene_4mer.xyz
set -euo pipefail

if [[ $# -lt 6 ]]; then
  sed -n '2,16p' "$0" >&2
  exit 2
fi
name="$1"; density="$2"; molar_mass="$3"; out_dir="$4"; shift 4

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bin="${BUILD_DIR:-$root/build}"
functional="${FUNCTIONAL:-CAM-B3LYP}"
basis="${BASIS:-def2-SVPD}"
aux="${AUX_BASIS:-def2-universal-JKFIT}"
grid="${GRID:-75,302}"
wavelengths="${WAVELENGTHS:-400,450,500,550,589,600,650,700,750,800,850,900,950,1000}"
gamma="${GAMMA:-0}"
states="${STATES:-0}"
order="${ORDER_PARAMETER:-0}"

mkdir -p "$out_dir"
results=()
for xyz in "$@"; do
  stem="$(basename "$xyz" .xyz)"
  result="$out_dir/$stem.json"
  if [[ ! -s "$result" ]]; then
    echo "response: $stem ($functional/$basis)" >&2
    aux_option=()
    if [[ -n "$aux" ]]; then aux_option=(--aux-basis "$aux"); fi
    "$bin/ellidft-response" "$xyz" "$basis" --functional "$functional" --grid "$grid" "${aux_option[@]}" \
      --wavelengths "$wavelengths" --gamma "$gamma" --states "$states" > "$result.tmp"
    mv "$result.tmp" "$result"
  else
    echo "response: $stem (kept existing $result)" >&2
  fi
  results+=("$result")
done

"$bin/ellidft-optics" --name "$name" --density "$density" --molar-mass "$molar_mass" \
  --order-parameter "$order" --table "$out_dir/$name.nk.txt" "${results[@]}" > "$out_dir/$name.optics.json"
echo "wrote $out_dir/$name.optics.json and $out_dir/$name.nk.txt" >&2
