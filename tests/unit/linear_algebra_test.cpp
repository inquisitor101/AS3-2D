#include "linear_algebra.hpp"

#include <gtest/gtest.h>

#include <initializer_list>
#include <limits>

namespace
{
  CMatrixAS3<as3double> Matrix(std::size_t rows, std::size_t columns,
                              std::initializer_list<as3double> values)
  {
    CMatrixAS3<as3double> result(rows, columns);
    std::size_t i = 0;
    for(const auto value : values) result[i++] = value;
    return result;
  }

  void ExpectMatrix(const CMatrixAS3<as3double> &actual,
                    const CMatrixAS3<as3double> &expected)
  {
    ASSERT_EQ(actual.row(), expected.row());
    ASSERT_EQ(actual.col(), expected.col());
    for(std::size_t i=0; i<actual.size(); i++)
      EXPECT_NEAR(actual[i], expected[i], 128.0*std::numeric_limits<as3double>::epsilon()) << i;
  }

  TEST(LinearAlgebra, RectangularMatrixProducts)
  {
    auto A = Matrix(2, 3, {1, 2, 3, 4, 5, 6});
    auto B = Matrix(3, 2, {7, 8, 9, 10, 11, 12});
    const auto expected = Matrix(2, 2, {58, 64, 139, 154});
    ExpectMatrix(NLinearAlgebra::MatrixMatrixMult(A, B), expected);
    ExpectMatrix(NLinearAlgebra::BLAS_GEMM(A, B), expected);
  }

  TEST(LinearAlgebra, MultipleVectorsStoredAsRows)
  {
    auto A = Matrix(2, 3, {1, 2, 3, 4, 5, 6});
    auto vectors = Matrix(2, 3, {1, -1, 2, 0, 3, -2});
    const auto expected = Matrix(2, 2, {5, 11, 0, 3});
    ExpectMatrix(NLinearAlgebra::MatrixVectorMult(A, vectors), expected);
    CMatrixAS3<as3double> actual(2, 2), temporary(2);
    NLinearAlgebra::MatrixVectorTransMult(A, vectors, actual, temporary);
    ExpectMatrix(actual, expected);
  }

  TEST(LinearAlgebra, LargeRectangularProducts)
  {
    // Preserve the original Skylabs 10x20 by 20x15 case with an analytical answer.
    constexpr std::size_t rows = 10, inner = 20, columns = 15;
    CMatrixAS3<as3double> A(rows, inner), B(inner, columns), expected(rows, columns);
    for(std::size_t i=0; i<rows; i++)
      for(std::size_t k=0; k<inner; k++) A(i,k) = static_cast<as3double>(i+k+1);
    for(std::size_t k=0; k<inner; k++)
      for(std::size_t j=0; j<columns; j++) B(k,j) = 2.0*static_cast<as3double>(k+1)-static_cast<as3double>(j);
    const as3double count = static_cast<as3double>(inner);
    const as3double sum = count*(count-1.0)/2.0;
    const as3double squares = count*(count-1.0)*(2.0*count-1.0)/6.0;
    for(std::size_t i=0; i<rows; i++)
      for(std::size_t j=0; j<columns; j++)
      {
        const as3double a = static_cast<as3double>(i+1), b = 2.0-static_cast<as3double>(j);
        expected(i,j) = count*a*b+(2.0*a+b)*sum+2.0*squares;
      }
    ExpectMatrix(NLinearAlgebra::MatrixMatrixMult(A, B), expected);
    ExpectMatrix(NLinearAlgebra::BLAS_GEMM(A, B), expected);
  }

  TEST(LinearAlgebra, AnalyticalThreeByThreeInverse)
  {
    auto A = Matrix(3, 3, {2, 1, 3, 0, 2, 4, 1, 1, 2});
    NLinearAlgebra::InverseMatrix(A);
    ExpectMatrix(A, Matrix(3, 3, {0, -0.5, 1, -2, -0.5, 4, 1, 0.5, -2}));
  }

  TEST(LinearAlgebra, TenByTenDiagonalInverse)
  {
    CMatrixAS3<as3double> A(10, 10), expected(10, 10);
    for(std::size_t i=0; i<10; i++)
    {
      A(i,i) = static_cast<as3double>(i+1);
      expected(i,i) = 1.0/A(i,i);
    }
    NLinearAlgebra::InverseMatrix(A);
    ExpectMatrix(A, expected);
  }

  TEST(LinearAlgebra, InverseWithPivoting)
  {
    auto original = Matrix(3, 3, {3, 0, 2, 2, 0, -2, 0, 1, 1});
    auto inverse = original;
    NLinearAlgebra::InverseMatrix(inverse);
    ExpectMatrix(inverse, Matrix(3, 3, {0.2, 0.2, 0, -0.2, 0.3, 1, 0.2, -0.3, 0}));
    const auto identity = NLinearAlgebra::CreateIdentityMatrix<as3double>(3);
    ExpectMatrix(NLinearAlgebra::MatrixMatrixMult(original, inverse), identity);
    ExpectMatrix(NLinearAlgebra::MatrixMatrixMult(inverse, original), identity);
  }

  TEST(LinearAlgebra, RectangularTranspose)
  {
    auto A = Matrix(2, 3, {1, 2, 3, 4, 5, 6});
    const auto expected = Matrix(3, 2, {1, 4, 2, 5, 3, 6});
    CMatrixAS3<as3double> result;
    NLinearAlgebra::TransposeAS3Matrix(A, result);
    ExpectMatrix(result, expected);
    ExpectMatrix(NLinearAlgebra::TransposeAS3Matrix(A), expected);
  }

  TEST(LinearAlgebra, KroneckerProduct)
  {
    auto A = Matrix(2, 2, {1, 2, 3, 4});
    auto B = Matrix(1, 2, {5, 6});
    ExpectMatrix(NLinearAlgebra::KroneckerProduct(A, B), Matrix(2, 4, {5, 6, 10, 12, 15, 18, 20, 24}));
  }

  TEST(LinearAlgebraDeathTest, RejectsSingularInverse)
  {
    auto singular = Matrix(2, 2, {1, 2, 2, 4});
    EXPECT_EXIT(NLinearAlgebra::InverseMatrix(singular), ::testing::ExitedWithCode(1), "");
  }

  TEST(LinearAlgebraDeathTest, RejectsIncompatibleProducts)
  {
    CMatrixAS3<as3double> A(2, 3), B(2, 2);
    EXPECT_EXIT((void)NLinearAlgebra::MatrixMatrixMult(A, B), ::testing::ExitedWithCode(1), "");
    EXPECT_EXIT((void)NLinearAlgebra::BLAS_GEMM(A, B), ::testing::ExitedWithCode(1), "");
  }
}
