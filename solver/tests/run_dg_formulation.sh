#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/as3-formulation-tests.XXXXXX")"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
omp_prefix="${OMP_PREFIX:-/opt/homebrew/opt/libomp}"
trap 'rm -rf "$build_dir"' EXIT

# Build the full solver with bounds checks, OpenMP and memory/undefined-behavior sanitizers.
if ! make -C "$project_dir" -f Makefile.clang MODE=debug OPENMP=1 CXX="$compiler" \
  OMP_PREFIX="$omp_prefix" BUILD_ROOT="$build_dir/build" BIN_DIR="$build_dir/bin" \
  EXTRA_CXXFLAGS='-fsanitize=address,undefined -fno-omit-frame-pointer' \
  EXTRA_LDFLAGS='-fsanitize=address,undefined' -j4 > "$build_dir/build.log" 2>&1; then
  cat "$build_dir/build.log"
  exit 1
fi
objects=()
while IFS= read -r source; do
  relative="${source#"$project_dir"/src/}"
  if [[ "$relative" != AS3.cpp ]]; then
    objects+=("$build_dir/build/debug-omp1-ltooff-native0/obj/${relative%.cpp}.o")
  fi
done < <(find -H "$project_dir/src" -type f -name '*.cpp' | sort)
"$compiler" -std=c++20 -O0 -g -DDEBUG -DENABLE_NAN_CHECK -DHAVE_OPENMP -fopenmp \
  -fsanitize=address,undefined -fno-omit-frame-pointer -Wall -Wextra -Wconversion \
  -I"$project_dir/include" -I"$omp_prefix/include" -L"$omp_prefix/lib" -Wl,-rpath,"$omp_prefix/lib" \
  "$project_dir/tests/dg_formulation_test.cpp" "${objects[@]}" -llapack -lblas \
  -o "$build_dir/dg_formulation_test"

# Absolute grid paths let the regression run without writing inside the example.
python3 - "$project_dir" "$build_dir" <<'PY'
from pathlib import Path
import re
import sys

project, target = map(Path, sys.argv[1:])
example = project.parent / 'example/curved_p5'
text = (example / 'param.cfg').read_text()
for key in ('TYPE_DG_METHOD', 'TYPE_DG_FORMULATION', 'TYPE_QUADRATURE', 'BONUS_QUADRATURE_FACTOR', 'GRID_DIRECTORY'):
    text = re.sub(rf'^[ \t]*{key}[ \t]*=.*\n?', '', text, flags=re.MULTILINE)
text += f'\nGRID_DIRECTORY = {example / "grid_lgl_p5"}\nBONUS_QUADRATURE_FACTOR = 0\n'
collocated = text + '\nTYPE_DG_METHOD = DG_SEM\nTYPE_QUADRATURE = LEGENDRE_GAUSS_LOBATTO\n'
(target / 'default.cfg').write_text(collocated)
(target / 'weak.cfg').write_text(collocated + '\nTYPE_DG_FORMULATION = WEAK_FORM\n')
(target / 'strong.cfg').write_text(collocated + '\nTYPE_DG_FORMULATION = STRONG_FORM\n')
(target / 'invalid.cfg').write_text(collocated + '\nTYPE_DG_FORMULATION = UNKNOWN_FORM\n')
for label, rule, bonus in [('GL', 'GAUSS_LEGENDRE', 0), ('LGL', 'LEGENDRE_GAUSS_LOBATTO', 0), ('OVER', 'GAUSS_LEGENDRE', 1)]:
    base = re.sub(r'^BONUS_QUADRATURE_FACTOR.*$', '', text, flags=re.MULTILINE)
    base += f'\nTYPE_DG_METHOD = NODAL_DG\nTYPE_QUADRATURE = {rule}\nBONUS_QUADRATURE_FACTOR = {bonus}\n'
    for form in ('WEAK_FORM', 'STRONG_FORM'):
        (target / f'{label}_{form}.cfg').write_text(base + f'\nTYPE_DG_FORMULATION = {form}\n')
PY

for threads in 1 4; do
  # Exercise both the absent option and explicit WEAK_FORM.
  weak_config="$build_dir/default.cfg"
  if [[ "$threads" == 4 ]]; then weak_config="$build_dir/weak.cfg"; fi
  if ! OMP_NUM_THREADS="$threads" OMP_DYNAMIC=FALSE "$build_dir/dg_formulation_test" \
    "$weak_config" "$build_dir/strong.cfg" > "$build_dir/run.log" 2>&1; then
    cat "$build_dir/run.log"
    exit 1
  fi
  rg 'checks passed|Imported-grid free-stream maximum' "$build_dir/run.log"
done
for label in GL LGL OVER; do
  for threads in 1 4; do
    if ! OMP_NUM_THREADS="$threads" OMP_DYNAMIC=FALSE "$build_dir/dg_formulation_test" \
      "$build_dir/${label}_WEAK_FORM.cfg" "$build_dir/${label}_STRONG_FORM.cfg" > "$build_dir/run.log" 2>&1; then
      cat "$build_dir/run.log"
      exit 1
    fi
    rg 'checks passed|Imported-grid free-stream maximum' "$build_dir/run.log"
  done
done
if "$build_dir/dg_formulation_test" "$build_dir/default.cfg" "$build_dir/invalid.cfg" > "$build_dir/error.log" 2>&1; then
  printf '%s\n' 'Unexpected success for invalid formulation.' >&2
  exit 1
fi
rg -q 'TYPE_DG_FORMULATION' "$build_dir/error.log"
printf '%s\n' 'All DG formulation checks passed.'
