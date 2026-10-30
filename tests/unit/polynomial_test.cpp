#include "polynomial_structure.hpp"
#include "quadrature/legendre_gauss_lobatto.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

namespace
{
  class Polynomial : public ::testing::TestWithParam<int> {};

  TEST_P(Polynomial, LagrangeInterpolationAndDifferentiation)
  {
    const auto n = static_cast<std::size_t>(GetParam());
    CMatrixAS3<as3double> nodes(n), evaluation(5), V, dV, Ve, dVe, L, D, De;
    const auto x = NLegendreGaussLobatto::GetPoints(n);
    for(std::size_t i=0; i<n; i++) nodes[i] = x[i];
    evaluation[0] = -1.0;
    evaluation[1] = -0.67;
    evaluation[2] = 0.0;
    evaluation[3] = 0.23;
    evaluation[4] = 1.0;

    NPolynomialUtility::OrthonormalLegendreBasis1D(n, nodes, V);
    NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D(n, nodes, dV);
    NPolynomialUtility::OrthonormalLegendreBasis1D(n, evaluation, Ve);
    NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D(n, evaluation, dVe);
    NPolynomialUtility::LagrangeBasis1D(Ve, V, L);
    NPolynomialUtility::DerivativeLagrangeBasis1D(dV, V, D);
    NPolynomialUtility::DerivativeLagrangeBasis1D(dVe, V, De);
    const as3double tolerance = 1024.0*static_cast<as3double>(n*n)
                               *std::numeric_limits<as3double>::epsilon();

    // Cardinality at the interpolation nodes.
    CMatrixAS3<as3double> cardinal;
    NPolynomialUtility::LagrangeBasis1D(V, V, cardinal);
    for(std::size_t i=0; i<n; i++)
      for(std::size_t j=0; j<n; j++) EXPECT_NEAR(cardinal(i,j), i == j ? 1.0 : 0.0, tolerance);

    // Analytical monomials verify partition of unity, differentiation of a
    // constant, and every polynomial the nodal basis must reproduce.
    for(std::size_t degree=0; degree<n; degree++)
    {
      SCOPED_TRACE(degree);
      for(std::size_t i=0; i<evaluation.size(); i++)
      {
        as3double value = 0.0, derivative = 0.0;
        for(std::size_t j=0; j<n; j++)
        {
          const as3double nodal = std::pow(nodes[j], static_cast<int>(degree));
          value += L(i,j)*nodal;
          derivative += De(i,j)*nodal;
        }
        const as3double exactDerivative = degree == 0 ? 0.0 :
          static_cast<as3double>(degree)*std::pow(evaluation[i], static_cast<int>(degree-1));
        EXPECT_NEAR(value, std::pow(evaluation[i], static_cast<int>(degree)), tolerance);
        EXPECT_NEAR(derivative, exactDerivative, tolerance);
      }
      for(std::size_t i=0; i<n; i++)
      {
        as3double derivative = 0.0;
        for(std::size_t j=0; j<n; j++) derivative += D(i,j)*std::pow(nodes[j], static_cast<int>(degree));
        const as3double exact = degree == 0 ? 0.0 :
          static_cast<as3double>(degree)*std::pow(nodes[i], static_cast<int>(degree-1));
        EXPECT_NEAR(derivative, exact, tolerance);
      }
    }

    CMatrixAS3<as3double> minimum, maximum;
    NPolynomialUtility::DerivativeLagrangeBasisFace1D(D, minimum, maximum);
    for(std::size_t j=0; j<n; j++)
    {
      EXPECT_NEAR(minimum[j], D(0,j), tolerance);
      EXPECT_NEAR(maximum[j], D(n-1,j), tolerance);
    }
  }

  INSTANTIATE_TEST_SUITE_P(PointCounts, Polynomial, ::testing::Range(2, 21));
}
