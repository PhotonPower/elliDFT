#!/usr/bin/env bash
# Configures, builds and tests elliDFT with one command.
#
#   scripts/build_and_test.sh            Release build in ./build, short tests only (a few seconds)
#   scripts/build_and_test.sh --slow     also the tests labelled "slow": full acceptance sets, minutes to
#                                        hours, meant for a server rather than the development machine
#   BUILD_DIR=build-debug BUILD_TYPE=Debug scripts/build_and_test.sh
#   scripts/build_and_test.sh -R unit    further arguments are passed to ctest
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${BUILD_DIR:-$root/build}"
build_type="${BUILD_TYPE:-Release}"

label_filter=(-LE slow)
if [[ "${1:-}" == "--slow" ]]; then
  label_filter=()
  shift
fi

generator=()
if command -v ninja >/dev/null 2>&1; then
  generator=(-G Ninja)
fi

cmake -S "$root" -B "$build_dir" "${generator[@]}" \
  -DCMAKE_BUILD_TYPE="$build_type" -DELLIDFT_WERROR=ON
cmake --build "$build_dir" --parallel
ctest --test-dir "$build_dir" --output-on-failure "${label_filter[@]}" "$@"
