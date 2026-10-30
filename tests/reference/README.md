# Reference test sources

There is one current AS3 solver. All runnable checks are regular tests in
`tests/unit/` and `tests/regression/`.

This directory preserves the original Skylabs test implementations before their
API updates. They are reference copies, excluded from compilation to avoid
duplicating the active tests. Their original Makefiles assume the former layout.
The active test runner is `../run_tests.sh`, backed by CMake and GoogleTest.

The checks have been incorporated as follows:

| Reference sources | Active checks |
| --- | --- |
| `test_unit_tensor_product.cpp` and `test_tensor_product.cpp` | `unit/tensor_product_test.cpp`: volume values/derivatives, all four surface values/derivatives and residual source/r/s/combined contributions across every compiled factory dimension |
| `test_unit_linear_algebra.cpp` | `unit/linear_algebra_test.cpp`: matrix checks, 10x20-by-20x15 products, analytical inverse and 10x10 diagonal inverse |
| `test_unit_standard_element.cpp` | `unit/quadrature_test.cpp` and `unit/polynomial_test.cpp`: solution nodes, quadrature, tensor-product weights, interpolation and differentiation |
| Original random/dense helpers | Deterministic inputs and dense 2D reference contractions in the active suites |

The original nine meshes and parameter file are in `../data/`. All 24 original
files are accounted for in `../reference_manifest.json`, with checksums.
The solver's DG-SEM/NODAL-DG work-in-progress behavior is unchanged by this
reorganization.
