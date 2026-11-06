#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/as3-collocated-tests.XXXXXX")"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
trap 'rm -rf "$build_dir"' EXIT
sources=(NError quadrature/CIntegrationRule linear_algebra/NLinearAlgebra
         NPolynomialUtility quadrature/CGaussJacobiQuadrature element/CStandardElement TensorProducts)
for i in "${!sources[@]}"; do
  sources[$i]="$project_dir/src/${sources[$i]}.cpp"
done
"$compiler" -std=c++20 -O1 -g -DDEBUG -Wall -Wextra -Wconversion -fopenmp-simd \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$project_dir/include" "$project_dir/tests/collocated_tensor_test.cpp" \
  "${sources[@]}" -llapack -lblas -o "$build_dir/collocated_tensor_test"
"$build_dir/collocated_tensor_test"
for invalid in null noncollocated wrong_size no_min_endpoint no_max_endpoint; do
  if "$build_dir/collocated_tensor_test" "$invalid" >"$build_dir/error.txt" 2>&1; then
    printf 'Unexpected success for %s\n' "$invalid" >&2
    exit 1
  fi
  case "$invalid" in
    null) rg -q 'Missing standard element' "$build_dir/error.txt" ;;
    noncollocated) rg -q 'requires matching solution and integration nodes' "$build_dir/error.txt" ;;
    wrong_size) rg -q 'Mismatch in tensor K dimension' "$build_dir/error.txt" ;;
    no_min_endpoint|no_max_endpoint) rg -q 'requires solution nodes at both reference endpoints' "$build_dir/error.txt" ;;
  esac
done
printf '%s\n' 'All collocated tensor checks passed.'
