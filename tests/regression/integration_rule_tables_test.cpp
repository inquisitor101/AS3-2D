/*!
 * @file
 * @brief Check all stored quadrature rules and rejection of invalid data.
 */

#include "quadrature/gauss_legendre.hpp"
#include "quadrature/legendre_gauss_lobatto.hpp"

#include <array>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>

static_assert(std::is_same_v<decltype(NGaussLegendre::GetPoints(1)),
                             std::span<const as3double>>);
static_assert(std::is_same_v<decltype(NGaussLegendre::GetWeights(1)),
                             std::span<const as3double>>);
static_assert(std::is_same_v<decltype(NLegendreGaussLobatto::GetPoints(2)),
                             std::span<const as3double>>);
static_assert(std::is_same_v<decltype(NLegendreGaussLobatto::GetWeights(2)),
                             std::span<const as3double>>);

namespace
{
  void Require(bool condition, const char *message)
  {
    if( !condition ) throw std::runtime_error(message);
  }

  void CheckInvalidRules(void)
  {
    using NIntegrationRuleValidation::ValidateRule;
    constexpr auto GL  = ETypeQuadrature::GAUSS_LEGENDRE;
    constexpr auto LGL = ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO;

    const auto points  = NGaussLegendre::GetPoints(2);
    const auto weights = NGaussLegendre::GetWeights(2);
    Require(!ValidateRule(GL, {}, {}), "Empty data accepted.");
    Require(!ValidateRule(GL, points, weights.first(1)), "Mismatched sizes accepted.");
    Require(!ValidateRule(GL, NLegendreGaussLobatto::GetPoints(2), weights),
            "Endpoints accepted in a GL rule.");
    Require(!ValidateRule(LGL, points, weights), "Missing LGL endpoints accepted.");

    std::array<as3double, 2> changed = {points[1], points[0]};
    Require(!ValidateRule(GL, changed, weights), "Unordered points accepted.");
    changed = {points[0], points[1] + 0.01};
    Require(!ValidateRule(GL, changed, weights), "Asymmetric points accepted.");
    changed = {-1.1, 1.1};
    Require(!ValidateRule(GL, changed, weights), "Out-of-range points accepted.");

    changed = {0.0, 2.0};
    Require(!ValidateRule(GL, points, changed), "Nonpositive weight accepted.");
    changed = {1.1, 0.9};
    Require(!ValidateRule(GL, points, changed), "Asymmetric weights accepted.");
    changed = {0.5, 0.5};
    Require(!ValidateRule(GL, points, changed), "Incorrect weight sum accepted.");

    // This wrong rule passes symmetry, point sum zero and weight sum two.
    // The degree-two moment must distinguish it from Gauss-Legendre.
    changed = {-0.5, 0.5};
    Require(!ValidateRule(GL, changed, weights), "Incorrect symmetric rule accepted.");

    // Likewise, LGL endpoints and symmetry do not determine the weights.
    const std::array<as3double, 3> wrongWeights = {0.5, 1.0, 0.5};
    Require(!ValidateRule(LGL, NLegendreGaussLobatto::GetPoints(3), wrongWeights),
            "Incorrect symmetric LGL weights accepted.");

    changed = {points[0], std::numeric_limits<as3double>::quiet_NaN()};
    Require(!ValidateRule(GL, changed, weights), "NaN point accepted.");
    changed = {1.0, std::numeric_limits<as3double>::infinity()};
    Require(!ValidateRule(GL, points, changed), "Infinite weight accepted.");
  }
}

int main()
{
  try
  {
    Require(NGaussLegendre::Validate(), "Invalid Gauss-Legendre table.");
    Require(NLegendreGaussLobatto::Validate(), "Invalid Legendre-Gauss-Lobatto table.");
    CheckInvalidRules();
    std::cout << "All 39 GL/LGL tables passed symmetry, sums, bounds, ordering, "
                 "positivity and polynomial exactness checks.\n"
              << "Invalid-data rejection checks passed.\n";
  }
  catch(const std::exception &error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
