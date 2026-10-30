#!/usr/bin/env bash
set -euo pipefail

tests_dir="$(cd "$(dirname "$0")" && pwd -P)"
build_dir="${AS3_TEST_BUILD_DIR:-$tests_dir/build}"
jobs="${AS3_TEST_JOBS:-4}"
mkdir -p "$build_dir"
build_dir="$(cd "$build_dir" && pwd -P)"

configure=(-S "$tests_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Debug)
cache="$build_dir/CMakeCache.txt"
if [[ -f "$cache" ]]; then
  cached_source="$(sed -n 's/^CMAKE_HOME_DIRECTORY:INTERNAL=//p' "$cache")"
  cached_build="$(sed -n 's/^CMAKE_CACHEFILE_DIR:INTERNAL=//p' "$cache")"
  if [[ "$cached_source" != "$tests_dir" || "$cached_build" != "$build_dir" ]]; then
    printf '%s\n' "Refreshing a build cache copied from another location."
    configure=(--fresh "${configure[@]}")
  fi
fi

cmake "${configure[@]}" "$@"
cmake --build "$build_dir" --parallel "$jobs"
ctest --test-dir "$build_dir" --output-on-failure --parallel "$jobs"
