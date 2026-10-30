#pragma once

/*!
 * @file
 * @brief Constant-evaluable checks for reference GL and LGL quadrature data.
 */

#include <cstddef>
#include <limits>
#include <span>

#include "option_structure.hpp"

namespace NIntegrationRuleValidation
{
  /*! @brief Absolute value usable in C++20 constant evaluation. */
  [[nodiscard]] inline constexpr as3double Absolute(as3double value)
  {
    return value < 0.0 ? -value : value;
  }

  /*!
   * @brief Verify dimensions, bounds, ordering, symmetry and polynomial moments.
   *
   * @param[in] typeQuadrature GL or LGL on [-1,+1] with unit integration weight.
   * @param[in] points ascending reference points.
   * @param[in] weights corresponding quadrature weights.
   *
   * @return true if the stored scalar values pass every check.
   * @note The tolerance scales with scalar epsilon and the number of points to
   *   allow rounding during accumulation. These checks validate stored values;
   *   they do not certify digits discarded when the source literals are rounded.
   */
  [[nodiscard]] inline constexpr bool ValidateRule(ETypeQuadrature typeQuadrature,
                                                   std::span<const as3double> points,
                                                   std::span<const as3double> weights)
  {
    const std::size_t n = points.size();
    if( n == 0 || weights.size() != n ) return false;

    std::size_t degreeExact = 0;
    switch( typeQuadrature )
    {
      case(ETypeQuadrature::GAUSS_LEGENDRE):
      {
        if( !(points.front() > -1.0 && points.back() < 1.0) ) return false;
        degreeExact = 2*n - 1;
        break;
      }
      case(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO):
      {
        if( n < 2 || points.front() != -1.0 || points.back() != 1.0 ) return false;
        degreeExact = 2*n - 3;
        break;
      }
      default: return false;
    }

    // Conservative allowance for summing n rounded scalar values.
    const as3double tolerance = 64.0 * std::numeric_limits<as3double>::epsilon()
                                    * static_cast<as3double>(n);
    as3double pointSum  = 0.0;
    as3double weightSum = 0.0;
    for(std::size_t i=0; i<n; i++)
    {
      const as3double x = points[i];
      const as3double w = weights[i];

      // The comparisons also reject NaNs and infinities.
      if( !(x >= -1.0 && x <= 1.0) ) return false;
      if( !(w > 0.0 && w <= std::numeric_limits<as3double>::max()) ) return false;
      if( i > 0 && !(points[i-1] < x) ) return false;
      if( Absolute(x + points[n-1-i]) > tolerance ) return false;
      if( Absolute(w - weights[n-1-i]) > tolerance ) return false;

      pointSum  += x;
      weightSum += w;
    }

    if( Absolute(pointSum) > tolerance ) return false;
    if( Absolute(weightSum - 2.0) > tolerance ) return false;
    if( n%2 && Absolute(points[n/2]) > tolerance ) return false;

    // Integrate x^k: zero for odd k, 2/(k+1) for even k.
    // Testing the full exactness degree checks points and weights together;
    // symmetry and the two sums alone can also hold for an incorrect rule.
    for(std::size_t degree=0; degree<=degreeExact; degree++)
    {
      as3double integral = 0.0;
      for(std::size_t i=0; i<n; i++)
      {
        as3double power = 1.0;
        for(std::size_t k=0; k<degree; k++) power *= points[i];
        integral += weights[i] * power;
      }
      const as3double exact = degree%2 ? 0.0 : 2.0/static_cast<as3double>(degree + 1);
      if( Absolute(integral - exact) > tolerance ) return false;
    }
    return true;
  }
} // namespace NIntegrationRuleValidation
