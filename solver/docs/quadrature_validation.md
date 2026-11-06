# Quadrature table validation

`include/quadrature/NGaussLegendre.hpp` and
`include/quadrature/NLegendreGaussLobatto.hpp`, together with
`include/quadrature/NNewtonCotes.hpp`, store the reference rules
using `as3double`. Their `GetPoints(n)` and `GetWeights(n)` functions return
`std::span<const as3double>` views. The existing decimal literals and alignment
are preserved.

`CIntegrationRule` uses these tables directly for GL rules with 1 through 20
points and LGL rules with 2 through 20 points. Larger rules continue through
`CGaussJacobiQuadrature`, with a `WARNING` identifying the rule and point count
and stating that support is not guaranteed. The generator's existing limits
still apply. Both paths build the same tensor-product 2D weights.

Select equidistant quadrature with `TYPE_QUADRATURE = NEWTON_COTES`, or construct
`CIntegrationRule(ETypeQuadrature::NEWTON_COTES, degree)`. The second argument is
the requested polynomial degree, as for GL and LGL. Newton-Cotes uses the
midpoint rule for one point and the closed rules, including both endpoints, for
2 through 8 points. Its nodal distribution is `ETypeNodalDistribution::EQD`.
The weights were obtained by integrating the equidistant Lagrange basis with
exact rational arithmetic, then rounding the results to 25 significant digits.

For `n` points the polynomial degree of exactness is `n` when `n` is odd and
`n-1` when `n` is even. The class chooses the smallest sufficient rule:

| Requested degree | Points selected | Degree of exactness |
| --- | --- | --- |
| 0–1 | 1 | 1 |
| 2–3 | 3 | 3 |
| 4–5 | 5 | 5 |
| 6–7 | 7 | 7 |

The even-point rules are also available through `NNewtonCotes::GetPoints(n)`
and `GetWeights(n)`, but provide no additional polynomial exactness over the
preceding odd-point rule. All eight tables have strictly positive weights.
The first negative-weight closed rule has 9 points. Counts above 8 are rejected
with `ERROR`; this deliberately excludes isolated higher positive-weight rules
such as the 10-point rule. `ComputePointCount1D` and the constructor both issue
`ERROR` for requested degrees above 7, before allocation or point-count
arithmetic. There is no higher-order fallback for Newton-Cotes.

For `NODAL_DG`, the requested degree is `(2 + BONUS_QUADRATURE_FACTOR)*POLY_ORDER`;
the default bonus therefore allows solution degrees up to 3. Full solver runs
also require the resulting solution/integration point counts to be present in
the existing factory specialization list. `DG_SEM` still requires LGL quadrature.

This limit preserves positive integration weights, rather than guaranteeing
stability of a complete discretization. Composite Newton-Cotes rules retain
low polynomial exactness as the point count grows. Discrete least-squares
quadrature can provide positive weights and higher exactness with more points,
but uses a different construction and is outside this rule's scope. See
[Glaubitz, Discrete least squares quadrature rules on equidistant and arbitrary points](https://archive.mpim-bonn.mpg.de/1088/1/preprint_2018_23.pdf)
for the Newton-Cotes sign limit and the least-squares alternative.

`CStandardElement::ComputeLocationDOFs1D` also takes LGL solution nodes directly
from the shared tables, replacing its separate hard-coded constants. Larger
counts use `CIntegrationRule` and its warning-producing fallback. Equidistant
nodes continue to be generated with uniform spacing. A single node is placed
at zero for either distribution. The current solution-node types are LGL and
equidistant; GL is available as an integration quadrature.

Each header checks every stored rule with `static_assert(Validate())` when it is
included. A failed check prevents compilation. The same checks can be requested
explicitly with `NGaussLegendre::Validate()` or
`NLegendreGaussLobatto::Validate()` or `NNewtonCotes::Validate()`.

The shared checker in `include/quadrature/NIntegrationRuleValidation.hpp` verifies:

- Matching, nonempty point and weight arrays and the minimum rule size.
- Finite, strictly increasing points in the reference interval, with endpoints
  excluded for GL and included for LGL and closed Newton-Cotes. The one-point
  Newton-Cotes rule must use the midpoint; other Newton-Cotes rules must have
  uniform spacing and no more than 8 points.
- Positive, finite weights.
- Symmetric point pairs and matching weight pairs, with a central zero for odd
  point counts.
- Point sum zero and weight sum two, within rounding tolerance.
- All polynomial moments through degree `2*n-1` for GL, `2*n-3` for LGL and
  `n-1+n%2` for Newton-Cotes:
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

Run the Newton-Cotes checks from the repository root:

```sh
bash solver/tests/run_newton_cotes.sh
```

These exercise all 8 stored rules, degree selection, 1D and tensor-product
polynomial moments, exactness metadata, standard-element construction and the
`ERROR` paths, including requests near the maximum `std::size_t` value. The
test runs with address and undefined-behavior sanitizers and floating-point
exception checks.

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
