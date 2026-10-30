#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")/../.." && pwd)"
solver_dir="$project_dir/solver"
build_dir="$(mktemp -d "${TMPDIR:-/tmp}/as3-element-tests.XXXXXX")"
compiler="${CXX:-/opt/homebrew/opt/llvm/bin/clang++}"
trap 'rm -rf "$build_dir"' EXIT

sources=(error_structure integration_rule_structure linear_algebra
         polynomial_structure quadrature_structure standard_element_structure)
for i in "${!sources[@]}"; do
	sources[$i]="$solver_dir/src/${sources[$i]}.cpp"
done

"$compiler" -std=c++20 -O1 -g -DDEBUG -DENABLE_NAN_CHECK \
	-I"$solver_dir/include" "$project_dir/tests/regression/element_dimensions_test.cpp" \
	"${sources[@]}" -llapack -lblas -o "$build_dir/element_dimensions_test"

test_binary="$build_dir/element_dimensions_test"
"$test_binary"

expect_success()
{
	local expected="$1"
	shift
	local output
	output="$("$test_binary" "$@" 2>&1)"
	if [[ "$output" != *"$expected"* ]]; then
		printf '%s\n' "$output" >&2
		exit 1
	fi
}

expect_failure()
{
	local expected="$1"
	shift
	local output
	if output="$("$test_binary" "$@" 2>&1)"; then
		printf 'Unexpected success: %s\n' "$*" >&2
		exit 1
	fi
	if [[ "$output" != *"$expected"* ]]; then
		printf '%s\n' "$output" >&2
		exit 1
	fi
}

# Exercise supported GL/LGL pairs, including the largest compiled GL rule.
expect_success "K=6 M=6" 5 0 GL NODAL_DG
expect_success "K=6 M=8" 5 1 GL NODAL_DG
expect_success "K=6 M=7" 5 0 LGL NODAL_DG
expect_success "K=10 M=14" 9 1 GL NODAL_DG

# Reference elements can use pairs and sizes absent from the factory catalogue.
expect_success "K=5 M=6" 4 0 LGL NODAL_DG
expect_success "K=6 M=9" 5 1 LGL NODAL_DG
expect_success "K=10 M=15" 9 1 LGL NODAL_DG
expect_success "K=11 M=11" 10 0 GL NODAL_DG
expect_success "K=2 M=17" 1 30 GL NODAL_DG

# Check reference-element restrictions independently of the run-time method guard.
expect_failure "DG-SEM requires zero bonus" 5 1 LGL DG_SEM
expect_failure "DG-SEM requires zero bonus" 5 0 GL DG_SEM
# A collocated reference element can be constructed; configured DG_SEM runs
# are rejected by CConfig on main.
expect_success "K=6 M=6" 5 0 LGL DG_SEM

printf '%s\n' "All element-dimension regression cases passed."
