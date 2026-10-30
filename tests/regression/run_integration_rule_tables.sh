#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/../.." && pwd)"
solver_dir="$project_dir/solver"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/as3-quadrature-tests.XXXXXX")"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
trap 'rm -rf "$build_dir"' EXIT

"$compiler" -std=c++20 -O1 -Wall -Wextra -Wconversion -Werror \
  -I"$solver_dir/include" "$project_dir/tests/regression/integration_rule_tables_test.cpp" \
  -o "$build_dir/integration_rule_tables_test"

"$build_dir/integration_rule_tables_test"
