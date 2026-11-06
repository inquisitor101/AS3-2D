#pragma once

/*!
 * @file
 * @brief Precomputed equidistant Newton-Cotes points and weights on [-1,+1].
 *
 * Counts range from 1 through 8; n is the number of points, not a degree.
 * The one-point rule is the midpoint rule. Other rules are closed and include
 * both endpoints. Weights are exact integrals of the Lagrange basis at the
 * equidistant points, computed with rational arithmetic and rounded to 25
 * significant digits, matching the GL/LGL reference tables.
 * An odd n integrates through degree n; an even n through degree n-1.
 * Every stored weight is positive. The 9-point rule has negative weights;
 * all counts above 8 are deliberately rejected, including isolated higher
 * counts with positive weights (such as 10), to retain a contiguous safe range.
 * Array entries use as3double. Unsuffixed literals are rounded to double
 * before conversion; as3double currently aliases double.
 * Points are strictly increasing; weights correspond to points by index.
 * Returned spans refer to immutable static arrays and allocate no memory.
 * Unsupported counts terminate with ERROR; there is no generated fallback.
 * Every stored rule is checked at compile time by static_assert(Validate()).
 */

#include <array>
#include <cstddef>
#include <span>

#include "common/AS3Options.hpp"
#include "quadrature/NIntegrationRuleValidation.hpp"

namespace NNewtonCotes
{
  inline constexpr std::size_t MinPoints = 1;
  inline constexpr std::size_t MaxPoints = 8;

  namespace detail
  {
    // 1-point Newton-Cotes rule (midpoint).
    inline constexpr std::array<as3double, 1> Points1 =
    {
       0.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 1> Weights1 =
    {
      2.000000000000000000000000e+0,
    };

    // 2-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 2> Points2 =
    {
      -1.000000000000000000000000e+0,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 2> Weights2 =
    {
      1.000000000000000000000000e+0,
      1.000000000000000000000000e+0,
    };

    // 3-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 3> Points3 =
    {
      -1.000000000000000000000000e+0,
       0.000000000000000000000000e+0,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 3> Weights3 =
    {
      3.333333333333333333333333e-1,
      1.333333333333333333333333e+0,
      3.333333333333333333333333e-1,
    };

    // 4-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 4> Points4 =
    {
      -1.000000000000000000000000e+0,
      -3.333333333333333333333333e-1,
       3.333333333333333333333333e-1,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 4> Weights4 =
    {
      2.500000000000000000000000e-1,
      7.500000000000000000000000e-1,
      7.500000000000000000000000e-1,
      2.500000000000000000000000e-1,
    };

    // 5-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 5> Points5 =
    {
      -1.000000000000000000000000e+0,
      -5.000000000000000000000000e-1,
       0.000000000000000000000000e+0,
       5.000000000000000000000000e-1,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 5> Weights5 =
    {
      1.555555555555555555555556e-1,
      7.111111111111111111111111e-1,
      2.666666666666666666666667e-1,
      7.111111111111111111111111e-1,
      1.555555555555555555555556e-1,
    };

    // 6-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 6> Points6 =
    {
      -1.000000000000000000000000e+0,
      -6.000000000000000000000000e-1,
      -2.000000000000000000000000e-1,
       2.000000000000000000000000e-1,
       6.000000000000000000000000e-1,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 6> Weights6 =
    {
      1.319444444444444444444444e-1,
      5.208333333333333333333333e-1,
      3.472222222222222222222222e-1,
      3.472222222222222222222222e-1,
      5.208333333333333333333333e-1,
      1.319444444444444444444444e-1,
    };

    // 7-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 7> Points7 =
    {
      -1.000000000000000000000000e+0,
      -6.666666666666666666666667e-1,
      -3.333333333333333333333333e-1,
       0.000000000000000000000000e+0,
       3.333333333333333333333333e-1,
       6.666666666666666666666667e-1,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 7> Weights7 =
    {
      9.761904761904761904761905e-2,
      5.142857142857142857142857e-1,
      6.428571428571428571428571e-2,
      6.476190476190476190476190e-1,
      6.428571428571428571428571e-2,
      5.142857142857142857142857e-1,
      9.761904761904761904761905e-2,
    };

    // 8-point Newton-Cotes rule (closed).
    inline constexpr std::array<as3double, 8> Points8 =
    {
      -1.000000000000000000000000e+0,
      -7.142857142857142857142857e-1,
      -4.285714285714285714285714e-1,
      -1.428571428571428571428571e-1,
       1.428571428571428571428571e-1,
       4.285714285714285714285714e-1,
       7.142857142857142857142857e-1,
       1.000000000000000000000000e+0,
    };

    inline constexpr std::array<as3double, 8> Weights8 =
    {
      8.692129629629629629629630e-2,
      4.140046296296296296296296e-1,
      1.531250000000000000000000e-1,
      3.459490740740740740740741e-1,
      3.459490740740740740740741e-1,
      1.531250000000000000000000e-1,
      4.140046296296296296296296e-1,
      8.692129629629629629629630e-2,
    };

  } // namespace detail

  /*! @brief Whether the requested point count is in the positive-weight range. */
  [[nodiscard]] inline constexpr bool IsAvailable(std::size_t n) noexcept
  {
    return n >= MinPoints && n <= MaxPoints;
  }

  /*! @brief Ascending points for an n-point rule on [-1,+1]. */
  [[nodiscard]] inline constexpr std::span<const as3double> GetPoints(std::size_t n)
  {
    switch (n)
    {
      case 1: return detail::Points1;
      case 2: return detail::Points2;
      case 3: return detail::Points3;
      case 4: return detail::Points4;
      case 5: return detail::Points5;
      case 6: return detail::Points6;
      case 7: return detail::Points7;
      case 8: return detail::Points8;
      default:
        ERROR("Newton-Cotes tables support 1 through 8 points with positive weights; "
              "the 9-point rule has negative weights and higher counts are not supported.");
    }
  }

  /*! @brief Weights corresponding to GetPoints(n). */
  [[nodiscard]] inline constexpr std::span<const as3double> GetWeights(std::size_t n)
  {
    switch (n)
    {
      case 1: return detail::Weights1;
      case 2: return detail::Weights2;
      case 3: return detail::Weights3;
      case 4: return detail::Weights4;
      case 5: return detail::Weights5;
      case 6: return detail::Weights6;
      case 7: return detail::Weights7;
      case 8: return detail::Weights8;
      default:
        ERROR("Newton-Cotes tables support 1 through 8 points with positive weights; "
              "the 9-point rule has negative weights and higher counts are not supported.");
    }
  }

  /*! @brief Verify every stored rule, including spacing and polynomial exactness. */
  [[nodiscard]] inline constexpr bool Validate(void)
  {
    for(std::size_t n=MinPoints; n<=MaxPoints; n++)
    {
      if( !NIntegrationRuleValidation::ValidateRule(ETypeQuadrature::NEWTON_COTES,
                                                    GetPoints(n), GetWeights(n)) )
        return false;
    }
    return true;
  }

  static_assert(Validate(), "Invalid precomputed Newton-Cotes points or weights.");

} // namespace NNewtonCotes
