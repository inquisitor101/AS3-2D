#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/.." && pwd)"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
build_dir="$project_dir/build/release-omp0-ltooff-native0"

make -C "$project_dir" -f Makefile.clang OPENMP=0 CXX="$compiler" -j4
objects=()
for object in "$build_dir"/obj/*.o; do
  if [[ "$(basename "$object")" != AS3.o ]]; then
    objects+=("$object")
  fi
done
"$compiler" -std=c++20 -O1 -g -Wall -Wextra -Wconversion \
  -I"$project_dir/include" "$project_dir/tests/tensor_factory_test.cpp" \
  "${objects[@]}" -llapack -lblas -o "$build_dir/tensor_factory_test"
"$build_dir/tensor_factory_test"
if "$build_dir/tensor_factory_test" invalid > "$build_dir/tensor_factory_error.log" 2>&1; then
  printf '%s\n' 'Unexpected success for noncollocated DG-SEM.' >&2
  exit 1
fi
rg -q 'DG-SEM tensor products require collocated nodes' "$build_dir/tensor_factory_error.log"

if "$build_dir/tensor_factory_test" gl_solution > "$build_dir/tensor_factory_gl_error.log" 2>&1; then
  printf '%s\n' 'Unexpected success for unsupported GL solution nodes.' >&2
  exit 1
fi
rg -q 'Gauss-Legendre solution nodes are not supported by the current face kernels' "$build_dir/tensor_factory_gl_error.log"

# Exercise both configuration options before standard-element and factory construction.
rg -v '^[[:space:]]*(TYPE_DG_METHOD|TYPE_QUADRATURE|BONUS_QUADRATURE_FACTOR)[[:space:]]*=' \
  "$project_dir/../example/curved_p5/param.cfg" > "$build_dir/tensor_factory_base.cfg"
for method in NODAL_DG DG_SEM; do
  config_file="$build_dir/tensor_factory_$method.cfg"
  cp "$build_dir/tensor_factory_base.cfg" "$config_file"
  printf '\nTYPE_DG_METHOD = %s\nBONUS_QUADRATURE_FACTOR = 0\n' "$method" >> "$config_file"
  if [[ "$method" == DG_SEM ]]; then
    printf '%s\n' 'TYPE_QUADRATURE = LEGENDRE_GAUSS_LOBATTO' >> "$config_file"
  else
    printf '%s\n' 'TYPE_QUADRATURE = GAUSS_LEGENDRE' >> "$config_file"
  fi
  "$build_dir/tensor_factory_test" config "$config_file"
done
