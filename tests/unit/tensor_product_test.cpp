/*!
 * @file
 * @brief Ports the Skylabs volume, surface and residual tensor checks to the current API.
 *
 * Dense two-dimensional contractions provide independent references for the
 * sum-factorized kernels. Inputs are deterministic and include all variables.
 */
#include "tensor_structure.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

// GoogleTest uses this to give each catalogue entry a readable CTest name.
void PrintTo(const CSupportedElementDimensions::CEntry &dimensions, std::ostream *output)
{
  *output << "K" << dimensions.mNSol1D << "M" << dimensions.mNInt1D << "N" << dimensions.mNVar;
}

namespace
{
  using Field = std::vector<as3double>;
  using Fields = std::array<Field, 3>;
  enum class Operation { Volume, Surface, ResidualVolume };
  // RuntimeVolume/RuntimeSurface currently have a 10-by-10 workspace limit.
  // Larger catalogue entries are checked through both compiled interfaces.

  Field Input(std::size_t size, as3double phase)
  {
    Field result(size);
    for(std::size_t i=0; i<size; i++)
      result[i] = std::sin(0.13*static_cast<as3double>(i+1)+phase)
                + std::cos(0.07*static_cast<as3double>(i+1)-phase);
    return result;
  }

  Fields Storage(std::size_t size, as3double value = 0.0)
  {
    return {Field(size, value), Field(size, value), Field(size, value)};
  }

  void ExpectField(const Field &actual, const Field &expected,
                   std::size_t count, std::size_t contractionSize)
  {
    ASSERT_GE(actual.size(), count);
    ASSERT_GE(expected.size(), count);
    for(std::size_t i=0; i<count; i++)
    {
      const as3double tolerance = 256.0*std::numeric_limits<as3double>::epsilon()
        *static_cast<as3double>(contractionSize)*std::max(as3double(1.0), std::abs(expected[i]));
      EXPECT_NEAR(actual[i], expected[i], tolerance) << "entry=" << i;
    }
    // Partially compiled and runtime kernels must respect their variable count.
    for(std::size_t i=count; i<actual.size(); i++) EXPECT_EQ(actual[i], 47.0) << i;
  }

  Fields DenseVolume(const CStandardElement &element, const Field &input, std::size_t variables)
  {
    const std::size_t K = element.GetnSol1D(), M = element.GetnInt1D();
    const auto &L = element.GetLagrangeInt1D();
    const auto &D = element.GetDerLagrangeInt1D();
    auto expected = Storage(variables*M*M);
    for(std::size_t v=0; v<variables; v++)
      for(std::size_t b=0; b<M; b++)
        for(std::size_t a=0; a<M; a++)
          for(std::size_t j=0; j<K; j++)
            for(std::size_t i=0; i<K; i++)
            {
              const auto value = input[v*K*K+j*K+i];
              const auto out = v*M*M+b*M+a;
              expected[0][out] += L(b,j)*L(a,i)*value;
              expected[1][out] += L(b,j)*D(a,i)*value;
              expected[2][out] += D(b,j)*L(a,i)*value;
            }
    return expected;
  }

  Fields DenseSurface(const CStandardElement &element, const Field &input,
                       std::size_t variables, int face)
  {
    const std::size_t K = element.GetnSol1D(), M = element.GetnInt1D();
    const auto &L = element.GetLagrangeInt1D();
    const auto &D = element.GetDerLagrangeInt1D();
    const auto &normalD = face%2 ? element.GetDerLagrangeMaxFace1D() : element.GetDerLagrangeMinFace1D();
    const std::size_t boundary = face%2 ? K-1 : 0;
    auto expected = Storage(variables*M);
    for(std::size_t v=0; v<variables; v++)
      for(std::size_t q=0; q<M; q++)
        for(std::size_t j=0; j<K; j++)
          for(std::size_t i=0; i<K; i++)
          {
            const auto value = input[v*K*K+j*K+i];
            const auto out = v*M+q;
            if(face < 2)
            {
              const as3double select = i == boundary ? 1.0 : 0.0;
              expected[0][out] += L(q,j)*select*value;
              expected[1][out] += L(q,j)*normalD[i]*value;
              expected[2][out] += D(q,j)*select*value;
            }
            else
            {
              const as3double select = j == boundary ? 1.0 : 0.0;
              expected[0][out] += L(q,i)*select*value;
              expected[1][out] += D(q,i)*select*value;
              expected[2][out] += L(q,i)*normalD[j]*value;
            }
          }
    return expected;
  }

  template<std::size_t K, std::size_t M, std::size_t N>
  void CheckVolume(CTensorProduct<K,M,N> &tensor, const CStandardElement &element)
  {
    const auto input = Input(N*K*K, 0.2);
    const auto expected = DenseVolume(element, input, N);
    for(unsigned int mask=0; mask<8; mask++)
    {
      SCOPED_TRACE(mask);
      for(int mode=0; mode<3; mode++)
      {
        if(mode == 2 && (K > 10 || M > 10)) continue;
        SCOPED_TRACE(mode);
        auto output = Storage(N*M*M, 47.0);
        auto *value = mask&1u ? output[0].data() : nullptr;
        auto *dr    = mask&2u ? output[1].data() : nullptr;
        auto *ds    = mask&4u ? output[2].data() : nullptr;
        const std::size_t variables = mode == 0 ? N : N-1;
        if(mode == 0) tensor.CompileTimeVolume(input.data(), value, dr, ds);
        else if(mode == 1) tensor.PartialCompileTimeVolume(variables, input.data(), value, dr, ds);
        else tensor.RuntimeVolume(K, variables, M, element.GetLagrangeInt1DTrans().data(),
                                  element.GetDerLagrangeInt1DTrans().data(), input.data(), value, dr, ds);
        for(unsigned int component=0; component<3; component++)
          if(mask & (1u<<component)) ExpectField(output[component], expected[component], variables*M*M, K*K);
      }
    }
  }

  template<std::size_t K, std::size_t M, std::size_t N>
  void CheckSurface(CTensorProduct<K,M,N> &tensor, const CStandardElement &element)
  {
    const auto input = Input(N*K*K, -0.3);
    using Fixed = void (ITensorProduct::*)(const as3double*, as3double*, as3double*, as3double*);
    using Partial = void (ITensorProduct::*)(std::size_t, const as3double*, as3double*, as3double*, as3double*);
    using Runtime = void (ITensorProduct::*)(std::size_t, std::size_t, std::size_t,
      const as3double*, const as3double*, const as3double*, const as3double*,
      as3double*, as3double*, as3double*);
    const std::array<Fixed,4> fixed = {&ITensorProduct::CompileTimeSurfaceIMIN, &ITensorProduct::CompileTimeSurfaceIMAX,
                                     &ITensorProduct::CompileTimeSurfaceJMIN, &ITensorProduct::CompileTimeSurfaceJMAX};
    const std::array<Partial,4> partial = {&ITensorProduct::PartialCompileTimeSurfaceIMIN, &ITensorProduct::PartialCompileTimeSurfaceIMAX,
                                         &ITensorProduct::PartialCompileTimeSurfaceJMIN, &ITensorProduct::PartialCompileTimeSurfaceJMAX};
    const std::array<Runtime,4> runtime = {&ITensorProduct::RuntimeSurfaceIMIN, &ITensorProduct::RuntimeSurfaceIMAX,
                                         &ITensorProduct::RuntimeSurfaceJMIN, &ITensorProduct::RuntimeSurfaceJMAX};
    for(std::size_t face=0; face<4; face++)
    {
      SCOPED_TRACE(face);
      const auto expected = DenseSurface(element, input, N, static_cast<int>(face));
      for(unsigned int mask=0; mask<8; mask++)
      {
        SCOPED_TRACE(mask);
        for(int mode=0; mode<3; mode++)
        {
          if(mode == 2 && (K > 10 || M > 10)) continue;
          SCOPED_TRACE(mode);
          auto output = Storage(N*M, 47.0);
          auto *value = mask&1u ? output[0].data() : nullptr;
          auto *dr    = mask&2u ? output[1].data() : nullptr;
          auto *ds    = mask&4u ? output[2].data() : nullptr;
          const std::size_t variables = mode == 0 ? N : N-1;
          if(mode == 0) (tensor.*fixed[face])(input.data(), value, dr, ds);
          else if(mode == 1) (tensor.*partial[face])(variables, input.data(), value, dr, ds);
          else
          {
            const auto &normalD = face%2 ? element.GetDerLagrangeMaxFace1D() : element.GetDerLagrangeMinFace1D();
            (tensor.*runtime[face])(K, variables, M, element.GetLagrangeInt1DTrans().data(),
              element.GetDerLagrangeInt1DTrans().data(), normalD.data(), input.data(), value, dr, ds);
          }
          for(unsigned int component=0; component<3; component++)
            if(mask & (1u<<component)) ExpectField(output[component], expected[component], variables*M, K*K);
        }
      }
    }
  }

  template<std::size_t K, std::size_t M, std::size_t N>
  void CheckResidualVolume(CTensorProduct<K,M,N> &tensor, const CStandardElement &element)
  {
    const Fields input = {Input(N*M*M, 0.1), Input(N*M*M, 0.7), Input(N*M*M, -0.4)};
    const auto &L = element.GetLagrangeInt1D();
    const auto &D = element.GetDerLagrangeInt1D();
    for(unsigned int mask=0; mask<8; mask++)
    {
      SCOPED_TRACE(mask);
      Field expected(N*K*K, 0.0), output(N*K*K, 47.0);
      for(std::size_t v=0; v<N; v++)
        for(std::size_t j=0; j<K; j++)
          for(std::size_t i=0; i<K; i++)
            for(std::size_t b=0; b<M; b++)
              for(std::size_t a=0; a<M; a++)
              {
                const auto in = v*M*M+b*M+a;
                const auto out = v*K*K+j*K+i;
                if(mask&1u) expected[out] += L(b,j)*L(a,i)*input[0][in];
                if(mask&2u) expected[out] += L(b,j)*D(a,i)*input[1][in];
                if(mask&4u) expected[out] += D(b,j)*L(a,i)*input[2][in];
              }
      tensor.CompileTimeResidualVolume(mask&1u ? input[0].data() : nullptr,
        mask&2u ? input[1].data() : nullptr, mask&4u ? input[2].data() : nullptr, output.data());
      ExpectField(output, expected, N*K*K, M*M);
    }
  }

  template<std::size_t K, std::size_t M, std::size_t N>
  void CheckCase(Operation operation)
  {
    // Match the existing factory catalogue without relying on the old config API.
    constexpr bool lobatto = K == 6 && M == 7;
    constexpr std::size_t bonus = M == K || lobatto ? 0 : (K == 2 ? 2 : 1);
    CStandardElement element(ETypeDOF::LGL, K-1, bonus,
      lobatto ? ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO : ETypeQuadrature::GAUSS_LEGENDRE,
      ETypeDGMethod::NODAL_DG);
    ASSERT_EQ(element.GetnInt1D(), M);
    CTensorProduct<K,M,N> tensor(&element);
    switch(operation)
    {
      case Operation::Volume: CheckVolume(tensor, element); break;
      case Operation::Surface: CheckSurface(tensor, element); break;
      case Operation::ResidualVolume: CheckResidualVolume(tensor, element); break;
    }
  }

  class TensorProduct : public ::testing::TestWithParam<CSupportedElementDimensions::CEntry>
  {
    protected:
      void Check(Operation operation)
      {
        const auto dimensions = GetParam();
        SCOPED_TRACE(::testing::Message() << "K=" << dimensions.mNSol1D << " M=" << dimensions.mNInt1D);
#define AS3_RUN_TENSOR_TEST(K,M,N) \
        if(dimensions.mNSol1D == K && dimensions.mNInt1D == M && dimensions.mNVar == N) \
        { CheckCase<K,M,N>(operation); return; }
        AS3_SUPPORTED_ELEMENT_SPECIALIZATIONS(AS3_RUN_TENSOR_TEST)
#undef AS3_RUN_TENSOR_TEST
        FAIL() << "Missing tensor test specialization.";
      }
  };

  TEST_P(TensorProduct, VolumeAndDerivatives) { Check(Operation::Volume); }
  TEST_P(TensorProduct, FourFacesAndDerivatives) { Check(Operation::Surface); }
  TEST_P(TensorProduct, ResidualVolume) { Check(Operation::ResidualVolume); }
  INSTANTIATE_TEST_SUITE_P(FactoryDimensions, TensorProduct,
    ::testing::ValuesIn(CSupportedElementDimensions::mCombinations));
}
