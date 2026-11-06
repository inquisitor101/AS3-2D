#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/.." && pwd)"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/as3-newton-cotes-tests.XXXXXX")"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
trap 'rm -rf "$build_dir"' EXIT

sources=(error_structure integration_rule_structure linear_algebra
         polynomial_structure quadrature_structure standard_element_structure)
for i in "${!sources[@]}"; do
  sources[$i]="$project_dir/src/${sources[$i]}.cpp"
done

"$compiler" -std=c++20 -O1 -g -DDEBUG -DENABLE_NAN_CHECK -Wall -Wextra -Wconversion \
  -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I"$project_dir/include" "$project_dir/tests/newton_cotes_test.cpp" \
  "${sources[@]}" -llapack -lblas -o "$build_dir/newton_cotes_test"

test_binary="$build_dir/newton_cotes_test"
"$test_binary"

expect_failure()
{
  local output
  if output="$("$test_binary" "$@" 2>&1)"; then
    printf 'Unexpected success: %s\n' "$*" >&2
    exit 1
  fi
  if [[ "$output" != *"9-point rule has negative weights"* ]]; then
    printf '%s\n' "$output" >&2
    exit 1
  fi
}

# Reject the first unsupported degree and very large degrees before arithmetic.
for mode in degree count; do
  for value in 8 9 10 100 max; do
    expect_failure "$mode" "$value"
  done
done

# Reject invalid table sizes, including the isolated positive 10-point rule.
for mode in points weights; do
  for value in 0 9 10 max; do
    expect_failure "$mode" "$value"
  done
done
expect_failure element 4

printf '%s\n' 'All Newton-Cotes regression cases passed.'
