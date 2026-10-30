#include "quadrature/gauss_legendre.hpp"
#include "quadrature/legendre_gauss_lobatto.hpp"
#include "integration_rule_structure.hpp"
#include "standard_element_structure.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <span>
#include <iostream>
#include <sstream>
#include <string>

namespace
{
  // NError::Warning writes to std::cout. Restore its buffer even on exceptions.
  class COutputCapture
  {
    std::ostringstream mOutput;
    std::streambuf *mOriginal;

    public:
      COutputCapture() : mOriginal(std::cout.rdbuf(mOutput.rdbuf())) {}
      ~COutputCapture() { std::cout.rdbuf(mOriginal); }
      COutputCapture(const COutputCapture&) = delete;
      COutputCapture& operator=(const COutputCapture&) = delete;
      std::string GetOutput() const { return mOutput.str(); }
  };

  void ExpectTabulatedData(const CIntegrationRule &rule,
                           std::span<const as3double> points,
                           std::span<const as3double> weights)
  {
    const std::size_t n = points.size();
    ASSERT_EQ(rule.GetnInt1D(), n);
    ASSERT_EQ(rule.GetnInt2D(), n*n);
    for(std::size_t i=0; i<n; i++)
    {
      // Exact equality verifies that construction uses the stored constants.
      EXPECT_EQ(rule.GetrInt1D()[i], points[i]) << i;
      EXPECT_EQ(rule.GetwInt1D()[i], weights[i]) << i;
      for(std::size_t j=0; j<n; j++)
        EXPECT_EQ(rule.GetwInt2D()[j*n+i], weights[i]*weights[j]);
    }
  }

  class Quadrature : public ::testing::TestWithParam<int> {};
  class LobattoQuadrature : public ::testing::TestWithParam<int> {};

  TEST_P(Quadrature, GaussLegendreMoments)
  {
    const auto n = static_cast<std::size_t>(GetParam());
    const auto x = NGaussLegendre::GetPoints(n);
    const auto w = NGaussLegendre::GetWeights(n);
    COutputCapture output;
    const CIntegrationRule rule(ETypeQuadrature::GAUSS_LEGENDRE, 2*(n-1));
    EXPECT_TRUE(output.GetOutput().empty());
    ExpectTabulatedData(rule, x, w);
    ASSERT_TRUE(NIntegrationRuleValidation::ValidateRule(ETypeQuadrature::GAUSS_LEGENDRE, x, w));
    for(std::size_t degree=0; degree<2*n; degree++)
    {
      SCOPED_TRACE(degree);
      as3double sum = 0.0;
      for(std::size_t i=0; i<n; i++) sum += w[i]*std::pow(x[i], static_cast<int>(degree));
      const as3double exact = degree%2 ? 0.0 : 2.0/static_cast<as3double>(degree+1);
      EXPECT_NEAR(sum, exact, 16.0*std::numeric_limits<as3double>::epsilon());
    }
  }

  TEST_P(LobattoQuadrature, LegendreGaussLobattoMoments)
  {
    const auto n = static_cast<std::size_t>(GetParam());
    const auto x = NLegendreGaussLobatto::GetPoints(n);
    const auto w = NLegendreGaussLobatto::GetWeights(n);
    COutputCapture output;
    const CIntegrationRule rule(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 2*(n-2));
    EXPECT_TRUE(output.GetOutput().empty());
    ExpectTabulatedData(rule, x, w);
    ASSERT_TRUE(NIntegrationRuleValidation::ValidateRule(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, x, w));

    // Solution nodes must use these exact constants as well.
    const CStandardElement element(ETypeDOF::LGL, n-1, 0,
                                    ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
    ASSERT_EQ(element.GetrSol1D().size(), n);
    for(std::size_t i=0; i<n; i++) EXPECT_EQ(element.GetrSol1D()[i], x[i]);
    EXPECT_TRUE(output.GetOutput().empty());
    for(std::size_t degree=0; degree<=2*n-3; degree++)
    {
      SCOPED_TRACE(degree);
      as3double sum = 0.0;
      for(std::size_t i=0; i<n; i++) sum += w[i]*std::pow(x[i], static_cast<int>(degree));
      const as3double exact = degree%2 ? 0.0 : 2.0/static_cast<as3double>(degree+1);
      EXPECT_NEAR(sum, exact, 16.0*std::numeric_limits<as3double>::epsilon());
    }
  }

  INSTANTIATE_TEST_SUITE_P(PointCounts, Quadrature, ::testing::Range(1, 21));
  INSTANTIATE_TEST_SUITE_P(PointCounts, LobattoQuadrature, ::testing::Range(2, 21));

  TEST(StandardElementNodes, SingleNodeIsTheMidpoint)
  {
    for(const auto distribution : {ETypeDOF::LGL, ETypeDOF::EQD})
    {
      const CStandardElement element(distribution, 0, 0,
                                      ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
      ASSERT_EQ(element.GetrSol1D().size(), 1u);
      EXPECT_EQ(element.GetrSol1D()[0], 0.0);
    }
  }

  TEST(StandardElementNodes, EquidistantNodesRemainUniform)
  {
    for(std::size_t n=2; n<=10; n++)
    {
      SCOPED_TRACE(n);
      const CStandardElement element(ETypeDOF::EQD, n-1, 0,
                                      ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
      const auto &points = element.GetrSol1D();
      ASSERT_EQ(points.size(), n);
      EXPECT_EQ(points[0], -1.0);
      EXPECT_EQ(points[n-1], 1.0);
      as3double sum = points[0];
      const as3double spacing = 2.0/static_cast<as3double>(n-1);
      for(std::size_t i=1; i<n; i++)
      {
        EXPECT_GT(points[i], points[i-1]);
        EXPECT_NEAR(points[i]-points[i-1], spacing, 8.0*std::numeric_limits<as3double>::epsilon());
        sum += points[i];
      }
      EXPECT_NEAR(sum, 0.0, 16.0*std::numeric_limits<as3double>::epsilon());
    }
  }

  TEST(StandardElementNodes, UntabulatedLobattoNodesUseTheFallback)
  {
    for(const std::size_t n : {21u, 22u})
    {
      SCOPED_TRACE(n);
      COutputCapture output;
      const CStandardElement element(ETypeDOF::LGL, n-1, 0,
                                      ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
      const auto warning = output.GetOutput();
      EXPECT_NE(warning.find("Legendre-Gauss-Lobatto quadrature with " + std::to_string(n) + " points"),
                std::string::npos);
      EXPECT_NE(warning.find("not guaranteed to be supported"), std::string::npos);

      const CIntegrationRule expected(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 2*(n-2));
      ASSERT_EQ(element.GetrSol1D().size(), n);
      for(std::size_t i=0; i<n; i++) EXPECT_EQ(element.GetrSol1D()[i], expected.GetrInt1D()[i]);
    }
  }

  TEST(QuadratureInputsDeathTest, UnsupportedTableSizes)
  {
    EXPECT_EXIT((void)NGaussLegendre::GetPoints(0), ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT((void)NGaussLegendre::GetWeights(21), ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT((void)NLegendreGaussLobatto::GetPoints(1), ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT((void)NLegendreGaussLobatto::GetWeights(21), ::testing::ExitedWithCode(1), "");
  }

  TEST(QuadratureConstruction, ComputedRulesBeyondTheTables)
  {
    for(const auto type : {ETypeQuadrature::GAUSS_LEGENDRE, ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO})
    {
      for(const std::size_t n : {21u, 32u})
      {
        SCOPED_TRACE(n);
        const bool gl = type == ETypeQuadrature::GAUSS_LEGENDRE;
        const std::size_t degree = 2*n - (gl ? 2u : 4u);
        COutputCapture output;
        const CIntegrationRule rule(type, degree);
        const std::string warning = output.GetOutput();
        EXPECT_NE(warning.find("Warning"), std::string::npos);
        EXPECT_NE(warning.find(gl ? "Gauss-Legendre" : "Legendre-Gauss-Lobatto"), std::string::npos);
        EXPECT_NE(warning.find(std::to_string(n) + " points"), std::string::npos);
        EXPECT_NE(warning.find("CGaussJacobiQuadrature"), std::string::npos);
        EXPECT_NE(warning.find("not guaranteed to be supported"), std::string::npos);
        EXPECT_EQ(warning.find("not tabulated"), warning.rfind("not tabulated"));
        ASSERT_EQ(rule.GetnInt1D(), n);
        EXPECT_TRUE(NIntegrationRuleValidation::ValidateRule(type,
          {rule.GetrInt1D().data(), rule.GetnInt1D()},
          {rule.GetwInt1D().data(), rule.GetnInt1D()}));
        const auto &volume = rule.GetwInt2D();
        as3double sum = 0.0;
        for(std::size_t i=0; i<volume.size(); i++) sum += volume[i];
        EXPECT_NEAR(sum, 4.0, 256.0*std::numeric_limits<as3double>::epsilon());
      }
    }
  }
}
