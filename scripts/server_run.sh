#!/usr/bin/env bash
# One command for the server: build, all tests, the run-time measurements the
# milestone reports still lack, and a report to send back.
#
#   scripts/server_run.sh                 build, short and slow tests, measurements on the styrene tetramer
#   scripts/server_run.sh --series        additionally the oligomer series of M6 (polystyrene 1-4, PMMA 1-4)
#   scripts/server_run.sh --with-direct   additionally the tetramer without density fitting (long)
#   scripts/server_run.sh --skip-slow     without the slow tests
#   scripts/server_run.sh --quick         plumbing test in seconds: water instead of the tetramer, no slow tests
#   OUT=results/server_1 scripts/server_run.sh
#
# Everything is written to one directory (default results/server_<date>);
# report.md there is the summary. A failing step does not stop the run; the
# report lists it.
set -uo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
quick=0; series=0; with_direct=0; skip_slow=0
for arg in "$@"; do
  case "$arg" in
    --quick) quick=1; skip_slow=1 ;;
    --series) series=1 ;;
    --with-direct) with_direct=1 ;;
    --skip-slow) skip_slow=1 ;;
    *) echo "unknown option $arg" >&2; sed -n '2,13p' "$0" >&2; exit 2 ;;
  esac
done
out="${OUT:-$root/results/server_$(date +%Y%m%d_%H%M%S)}"
build="${BUILD_DIR:-$root/build}"
mkdir -p "$out/bench"
report="$out/report.md"
cd "$root"

say() { echo "== $*" >&2; }

# --- system ----------------------------------------------------------------
say "system"
{
  echo "# elliDFT: Lauf auf dem Server"
  echo
  echo "| | |"
  echo "|---|---|"
  echo "| Datum | $(date '+%Y-%m-%d %H:%M:%S') |"
  echo "| Rechner | $(uname -srm) |"
  if command -v lscpu >/dev/null 2>&1; then
    echo "| Prozessor | $(lscpu | sed -n 's/^Model name: *//p' | head -1) |"
    echo "| Kerne | $(nproc) |"
    echo "| Speicher | $(free -g | awk '/^Mem:/ {print $2 " GB"}') |"
  else
    echo "| Prozessor | $(sysctl -n machdep.cpu.brand_string 2>/dev/null) |"
    echo "| Kerne | $(sysctl -n hw.ncpu 2>/dev/null) |"
    echo "| Speicher | $(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 )) GB |"
  fi
  echo "| Übersetzer | $(${CXX:-c++} --version 2>/dev/null | head -1) |"
  echo "| CMake | $(cmake --version | head -1) |"
  echo "| Commit | $(git rev-parse --short HEAD 2>/dev/null) ($(git status --short 2>/dev/null | wc -l | tr -d ' ') geänderte Dateien) |"
  echo "| OMP_NUM_THREADS | ${OMP_NUM_THREADS:-nicht gesetzt} |"
  if command -v nvidia-smi >/dev/null 2>&1; then
    echo "| Grafikkarte | $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | head -1) |"
  else
    echo "| Grafikkarte | keine NVIDIA-Karte gefunden |"
  fi
  echo
} > "$report"

# --- build -------------------------------------------------------------------
say "build"
generator=()
command -v ninja >/dev/null 2>&1 && generator=(-G Ninja)
build_ok=ja
{ cmake -S "$root" -B "$build" "${generator[@]}" -DCMAKE_BUILD_TYPE=Release -DELLIDFT_WERROR=ON &&
  cmake --build "$build" --parallel; } > "$out/build.log" 2>&1 || build_ok=NEIN
{
  echo "## Übersetzen"
  echo
  echo "Erfolgreich: $build_ok (Protokoll: \`build.log\`)"
  echo
} >> "$report"
if [[ "$build_ok" != ja ]]; then
  say "build failed, see $out/build.log"
  exit 1
fi

# --- tests -------------------------------------------------------------------
run_tests() {   # name, ctest arguments ...
  local name="$1"; shift
  say "tests: $name"
  local start=$SECONDS
  ctest --test-dir "$build" --output-on-failure --output-junit "$out/tests_$name.xml" "$@" > "$out/tests_$name.log" 2>&1
  local status=$?
  local summary
  summary="$(grep -E "tests passed|tests failed" "$out/tests_$name.log" | tail -1)"
  echo "| $name | ${summary:-kein Ergebnis} | $((SECONDS - start)) s | $([[ $status -eq 0 ]] && echo ja || echo NEIN) |" >> "$report"
  if [[ $status -ne 0 ]]; then
    grep -E "\*\*\*Failed|\(Failed\)" "$out/tests_$name.log" | sed 's/^/    /' > "$out/tests_${name}_failed.txt"
  fi
}
{
  echo "## Tests"
  echo
  echo "| Lauf | Ergebnis | Dauer | bestanden |"
  echo "|---|---|---|---|"
} >> "$report"
run_tests kurz -LE slow
if [[ $skip_slow -eq 0 ]]; then
  run_tests slow -L slow
else
  echo "| slow | übersprungen | | |" >> "$report"
fi
{
  echo
  echo "Protokolle: \`tests_*.log\`, fehlgeschlagene Tests in \`tests_*_failed.txt\`."
  echo "Die Zeilen mit Messwerten der Vergleiche stehen in den Protokollen (suche nach \`relative deviation\`, \`|dE|\`)."
  echo
} >> "$report"

# --- measurements ------------------------------------------------------------
# Peak memory through GNU time (Linux) or BSD time (macOS), if available.
timed() {   # name, command ...
  local name="$1"; shift
  say "measurement: $name"
  local start; start=$(date +%s.%N 2>/dev/null || date +%s)
  local memory="-"
  if /usr/bin/time -f "%M" true >/dev/null 2>&1; then
    /usr/bin/time -f "%M" -o "$out/bench/$name.mem" "$@" > "$out/bench/$name.json" 2> "$out/bench/$name.log"
    local status=$?
    memory="$(awk '{printf "%.2f GB", $1 / 1048576}' "$out/bench/$name.mem" 2>/dev/null)"
  elif /usr/bin/time -l true >/dev/null 2>&1; then
    /usr/bin/time -l "$@" > "$out/bench/$name.json" 2> "$out/bench/$name.log"
    local status=$?
    memory="$(awk '/maximum resident set size/ {printf "%.2f GB", $1 / 1073741824}' "$out/bench/$name.log")"
  else
    "$@" > "$out/bench/$name.json" 2> "$out/bench/$name.log"
    local status=$?
  fi
  local stop; stop=$(date +%s.%N 2>/dev/null || date +%s)
  local wall; wall=$(awk -v a="$start" -v b="$stop" 'BEGIN {printf "%.1f", b - a}')
  local energy; energy="$(grep -m1 -E '"(energy_total|energy_hartree)"' "$out/bench/$name.json" | sed -E 's/.*: *(-?[0-9.]+([eE][-+]?[0-9]+)?).*/\1/')"
  local iterations; iterations="$(grep -m1 -E '"(iterations|scf_iterations)"' "$out/bench/$name.json" | sed -E 's/.*: *([0-9]+).*/\1/')"
  echo "| $name | $wall s | $memory | ${energy:--} | ${iterations:--} | $([[ $status -eq 0 ]] && echo ja || echo NEIN) |" >> "$report"
}

if [[ $quick -eq 1 ]]; then
  molecule=h2o; basis=STO-3G; jfit=def2-universal-jfit; jkfit=def2-universal-jkfit
else
  molecule=styrene_4mer; basis=def2-SVP; jfit=def2-universal-jfit; jkfit=def2-universal-jkfit
fi
geometry="$root/data/geometries/$molecule.xyz"
{
  echo "## Zeitmessungen"
  echo
  echo "Molekül \`$molecule\`, Basis $basis. Einzelheiten je Zeile in \`bench/<Name>.json\` und \`.log\`."
  echo
  echo "| Messung | Dauer | Speicherspitze | Energie / Hartree | Iterationen | gelaufen |"
  echo "|---|---|---|---|---|---|"
} >> "$report"
# M3: PBE with density-fitted Coulomb, against the reference program
timed "M3_PBE_RI-J" "$build/ellidft-scf" "$geometry" "$basis" --functional PBE --aux-basis "$jfit"
python="$root/.venv/bin/python"
if [[ -x "$python" ]] && "$python" -c "import pyscf" 2>/dev/null; then
  ( cd "$root/tests/reference" && timed "M3_PBE_RI-J_PySCF" "$python" time_scf.py "$molecule" "$basis" PBE --aux "$jfit" )
else
  echo "| M3_PBE_RI-J_PySCF | nicht gemessen: keine .venv mit PySCF (siehe docs/progress/README.md) | | | | |" >> "$report"
fi
# M4: hybrid functional, exchange with density fitting; direct exchange on request
timed "M4_PBE0_RI-JK" "$build/ellidft-scf" "$geometry" "$basis" --functional PBE0 --aux-basis "$jkfit"
if [[ $with_direct -eq 1 || $quick -eq 1 ]]; then
  timed "M3_PBE_direkt" "$build/ellidft-scf" "$geometry" "$basis" --functional PBE
  timed "M4_PBE0_direkt" "$build/ellidft-scf" "$geometry" "$basis" --functional PBE0
else
  echo "| M3_PBE_direkt, M4_PBE0_direkt | übersprungen (\`--with-direct\`) | | | | |" >> "$report"
fi
{
  echo
  echo "Kriterium M3: Dauer von \`M3_PBE_RI-J\` höchstens zehnmal die von \`M3_PBE_RI-J_PySCF\`."
  echo "M4: Vergleich von \`M4_PBE0_RI-JK\` mit \`M4_PBE0_direkt\` (Dauer, Speicher) für die Entscheidung über den Austauschaufbau."
  echo
} >> "$report"

# --- oligomer series (M6) ------------------------------------------------------
if [[ $series -eq 1 || $quick -eq 1 ]]; then
  {
    echo "## Oligomer-Reihen (M6)"
    echo
    echo "Dichten: Richtwerte des Auftrags (1,05 und 1,18 g/cm³), nicht die der Proben."
    echo
    echo "| Reihe | Dauer | gelaufen | Ergebnis |"
    echo "|---|---|---|---|"
  } >> "$report"
  run_series() {   # name density molar_mass geometries ...
    local name="$1"; shift
    say "series: $name"
    local start=$SECONDS
    BUILD_DIR="$build" "$root/scripts/run_oligomer_series.sh" "$name" "$@" > "$out/series_$name.log" 2>&1
    local status=$?
    echo "| $name | $((SECONDS - start)) s | $([[ $status -eq 0 ]] && echo ja || echo NEIN) | \`series_$name/$name.optics.json\`, \`.nk.txt\` |" >> "$report"
  }
  g="$root/data/geometries"
  if [[ $quick -eq 1 ]]; then
    FUNCTIONAL=HF BASIS=STO-3G AUX_BASIS=def2-universal-JKFIT WAVELENGTHS=589 \
      run_series quick 1.0 2.0 "$out/series_quick" "$g/h2.xyz" "$g/h2.xyz"
  else
    run_series polystyrol 1.05 104.15 "$out/series_polystyrol" "$g/ethylbenzene.xyz" "$g/styrene_2mer.xyz" "$g/styrene_3mer.xyz" "$g/styrene_4mer.xyz"
    run_series pmma 1.18 100.12 "$out/series_pmma" "$g/methyl_isobutyrate.xyz" "$g/mma_2mer.xyz" "$g/mma_3mer.xyz" "$g/mma_4mer.xyz"
  fi
  echo >> "$report"
fi

echo "Vollständig durchgelaufen: $(date '+%Y-%m-%d %H:%M:%S')" >> "$report"
say "done: $report"
