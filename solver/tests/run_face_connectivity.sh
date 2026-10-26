#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/as3-face-tests.XXXXXX")"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
omp_prefix="${OMP_PREFIX:-/opt/homebrew/opt/libomp}"

# Geometry and mesh import only; the full solver is tested separately.
sources=(config_structure error_structure face_structure geometry_structure
         import_structure input_structure log_structure marker_structure
         openmp_structure load_balance_faces param_structure option_structure)
for i in "${!sources[@]}"; do
  sources[$i]="$project_dir/src/${sources[$i]}.cpp"
done

"$compiler" -std=c++20 -O1 -g -DDEBUG -DHAVE_OPENMP -fopenmp \
  -I"$project_dir/include" -I"$omp_prefix/include" \
  -L"$omp_prefix/lib" -Wl,-rpath,"$omp_prefix/lib" -Wl,-dead_strip \
  "$project_dir/tests/face_connectivity_test.cpp" "${sources[@]}" \
  -o "$build_dir/face_connectivity_test"

# Keep mesh paths relative to the selected example, as in its param.cfg.
example_dir="${1:-$project_dir/../example/curved_p5}"
cd "$example_dir"
for threads in 1 4; do
  OMP_NUM_THREADS="$threads" "$build_dir/face_connectivity_test" param.cfg
done

# Exercise physical boundary collections and all three category offsets too.
sed '/^MARKER_BC_PERIODIC/d' param.cfg > "$build_dir/physical.cfg"
for edge in bottom top left right; do
  printf 'MARKER_BC_SYMMETRY = ( %s )\n' "$edge" >> "$build_dir/physical.cfg"
done
for threads in 1 4; do
  OMP_NUM_THREADS="$threads" "$build_dir/face_connectivity_test" "$build_dir/physical.cfg"
done

sed '/^MARKER_BC_PERIODIC = ( left/d' param.cfg > "$build_dir/mixed.cfg"
for edge in left right; do
  printf 'MARKER_BC_SYMMETRY = ( %s )\n' "$edge" >> "$build_dir/mixed.cfg"
done
for threads in 1 4; do
  OMP_NUM_THREADS="$threads" "$build_dir/face_connectivity_test" "$build_dir/mixed.cfg"
done
