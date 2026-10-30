# AS3 numerical tests

Tests live beside the solver. GoogleTest provides assertions, parameterized
cases and failure diagnostics; CMake builds the test executables and CTest runs
the discovered cases. The solver's GNU, Intel and Clang Makefiles remain separate.

```text
solver/
  include/quadrature/
    gauss_legendre.hpp
    legendre_gauss_lobatto.hpp
    integration_rule_validation.hpp
tests/
  CMakeLists.txt
  unit/
    quadrature_test.cpp
    polynomial_test.cpp
    linear_algebra_test.cpp
    tensor_product_test.cpp
  regression/
    ...existing standalone tests and shell runners...
  data/
    grid_AS3/
    param_t1.cfg
  reference/
    ...original test source copies...
  reference_manifest.json
```

The quadrature headers retain their compile-time `static_assert` validation and
`as3double` storage. Separate tests exercise the routines, corrupted inputs and
mathematical properties. Existing solver input checks and numerical checks remain
in place.

## Build and run

From the repository root, with CMake 3.24+, a C++20 compiler and LP64 BLAS/LAPACK:

```sh
cmake -S tests -B tests/build -DCMAKE_BUILD_TYPE=Debug
cmake --build tests/build -j4
ctest --test-dir tests/build --output-on-failure
```

On macOS, select Homebrew LLVM explicitly if desired by adding
`-DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++` to configuration.
CMake uses an installed GoogleTest package when available. Otherwise it downloads
the pinned GoogleTest 1.17.0 archive and verifies its SHA-256. For an offline
build, provide an installed package through `CMAKE_PREFIX_PATH`, or configure with
`-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=/absolute/path/to/googletest-1.17.0`.
Dependencies are test-only; no GoogleTest installation is required to build AS3.

Filter the discovered tests or run GoogleTest directly:

```sh
ctest --test-dir tests/build -R 'Quadrature|Lobatto' --output-on-failure
tests/build/as3_unit_tests --gtest_filter='LinearAlgebra*.*'
```

## Current coverage

| Area | Checks |
| --- | --- |
| Quadrature | Exact table use without warnings for every stored GL/LGL rule, full polynomial exactness, invalid counts and invalid data, computed 21- and 32-point rules with fallback warnings, tensor-product weights |
| Standard-element nodes | Exact shared LGL constants, generated 21- and 22-point LGL nodes, uniform equidistant spacing, single-node midpoint |
| Lagrange interpolation | Cardinality, partition of unity, reproduction of every monomial through the basis degree, at nodes and other evaluation points |
| Differentiation | Derivative of a constant and analytical monomial derivatives at nodes, endpoints and other evaluation points |
| Linear algebra | Known rectangular products for manual and BLAS routines, multiple row-stored vectors, inverse with pivoting and both inverse residuals, transpose, Kronecker product |
| Migrated Skylabs cases | Original larger matrix-product dimensions, analytical inverse and 10x10 diagonal inverse |
| Tensor products | Volume and all four surface values/derivatives, residual source/r/s/combined contributions, every factory dimension, compiled/partial/runtime interfaces and optional outputs |
| Errors | Singular matrix inversion and incompatible matrix products terminate as expected |
| Existing regressions | Quadrature table rejection checks and standard-element dimensions |

Expected answers use analytical polynomials and independently specified small
matrices. Avoid testing an implementation only against another routine that
shares its code. Use scalar-aware tolerances and include point count, degree or
matrix indices in diagnostics. Passing these numerical tests validates stored
precision, not all 25 digits printed in the quadrature source literals.

Add a suite under `unit/`, add its source to `as3_unit_tests` in `CMakeLists.txt`,
and define `TEST`, `TEST_F` or `TEST_P` cases. CTest discovers them automatically.
Additional mesh/solver integration cases can be added as separate targets without
linking the solver's `AS3.cpp` entry point.

## Current tests and reference material

All executable numerical checks live in `unit/` and `regression/`. They test the
current solver. The tensor and matrix checks brought over from Skylabs are
included in the same GoogleTest executable as the other unit tests.

Original source copies are retained in `reference/` to record the checks before
API adaptation. They are excluded from compilation to avoid maintaining duplicate
implementations of the same tests. Mesh fixtures and their original parameter
file live in `data/`. The original-to-current file mappings and checksums are in
`reference_manifest.json`. See [the coverage mapping](reference/README.md).

Run the suite from any directory with the portable wrapper:

```sh
bash /path/to/AS3-2D/tests/run_tests.sh
```

It derives source paths from its own location and uses CMake's normal compiler
selection (`CXX` or `-DCMAKE_CXX_COMPILER=...` can override it). If a copied build
cache points to another source/build location, it refreshes that cache before
configuration. Pass additional CMake options as arguments; `AS3_TEST_BUILD_DIR`
and `AS3_TEST_JOBS` control the build directory and parallelism. Generated build
artifacts still should not be transferred between checkouts or computers.

## Standalone regressions

The original runners still work from any working directory:

```sh
bash tests/regression/run_integration_rule_tables.sh
bash tests/regression/run_element_dimensions.sh
bash tests/regression/run_face_connectivity.sh /absolute/path/to/example/curved_p5
```

The face-connectivity runner additionally needs OpenMP/libomp and the example
mesh files. Its default fixture location is `example/curved_p5` beside `solver/`.
It is kept outside the default CTest run because fixtures are not part of this
workspace copy.

Framework references: [GoogleTest's CMake quickstart](https://google.github.io/googletest/quickstart-cmake.html)
and [CTest GoogleTest discovery](https://cmake.org/cmake/help/latest/module/GoogleTest.html).
