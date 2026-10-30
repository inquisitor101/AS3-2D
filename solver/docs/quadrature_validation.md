# Quadrature table validation

`include/quadrature/gauss_legendre.hpp` and
`include/quadrature/legendre_gauss_lobatto.hpp` store the reference rules
using `as3double`. Their `GetPoints(n)` and `GetWeights(n)` functions return
`std::span<const as3double>` views. The existing decimal literals and alignment
are preserved.

`CIntegrationRule` uses these tables directly for GL rules with 1 through 20
points and LGL rules with 2 through 20 points. Larger rules continue through
`CGaussJacobiQuadrature`, with a `WARNING` identifying the rule and point count
and stating that support is not guaranteed. The generator's existing limits
still apply. Both paths build the same tensor-product 2D weights.

`CStandardElement::ComputeLocationDOFs1D` also takes LGL solution nodes directly
from the shared tables, replacing its separate hard-coded constants. Larger
counts use `CIntegrationRule` and its warning-producing fallback. Equidistant
nodes continue to be generated with uniform spacing. A single node is placed
at zero for either distribution. The current solution-node types are LGL and
equidistant; GL is available as an integration quadrature.

Each header checks every stored rule with `static_assert(Validate())` when it is
included. A failed check prevents compilation. The same checks can be requested
explicitly with `NGaussLegendre::Validate()` or
`NLegendreGaussLobatto::Validate()`.

The shared checker in `include/quadrature/integration_rule_validation.hpp` verifies:

- Matching, nonempty point and weight arrays and the minimum rule size.
- Finite, strictly increasing points in the reference interval, with endpoints
  excluded for GL and included for LGL.
- Positive, finite weights.
- Symmetric point pairs and matching weight pairs, with a central zero for odd
  point counts.
- Point sum zero and weight sum two, within rounding tolerance.
- All polynomial moments through degree `2*n-1` for GL and `2*n-3` for LGL:
  the weighted sum of `x^k` must be zero for odd `k`, and `2/(k+1)` for even `k`.

Symmetry and the two sums alone do not establish that a rule is correct. For
example, points `{-0.5, 0.5}` with weights `{1, 1}` satisfy both sums but fail to
integrate `x^2`. The regression test checks rejection of this case, along with
invalid bounds, ordering, weights, dimensions and nonfinite values.

The absolute tolerance is `64*n*std::numeric_limits<as3double>::epsilon()` to
allow rounding in arithmetic and accumulation. These checks validate the stored
scalar values; they do not certify all 25 significant digits in the source text.
Unsuffixed floating literals are rounded to `double` before conversion to
`as3double`, which currently aliases `double`.

The standalone regressions live in `tests/` beside `solver/`. Run all 39 stored
rules and the invalid-data cases from the repository root:

```sh
bash tests/regression/run_integration_rule_tables.sh
```

The decimal text was independently compared against SymPy's GL/LGL routines at
60-digit precision, with all 838 values agreeing after rounding to 25 significant
digits. The GL constants for 2 through 20 points were also compared against GSL's
published tables and agreed at their printed precision and as doubles. See the
[SymPy quadrature documentation](https://docs.sympy.org/latest/modules/integrals/integrals.html#numeric-integrals)
and [GSL reference constants](https://raw.githubusercontent.com/ampl/gsl/master/integration/glfixed.c).

The GoogleTest suites also verify quadrature, Lagrange interpolation and
differentiation, and linear algebra. See [the test setup](../../tests/README.md)
for CMake/CTest commands and adding more suites.
