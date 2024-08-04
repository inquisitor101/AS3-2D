#pragma once

#include "option_structure.hpp"


/*!
 * @brief A namespace used for storing linear algebra functionalities.
 */
namespace NLinearAlgebra
{
	/*!
	 * @brief Function that computes the inverse of the given matrix (double-precision).
	 *
	 * @param[in] n dimension of the (square) matrix.
	 * @param[in,out] A inverse of the matrix. 
	 */
	void InverseMatrix(size_t nn, CMatrixAS3<as3double> &A);

	/*!
	 * @brief Function that computes a generalized matrix-matrix multiplication via BLAS (double-precision).
	 *
	 * @param[in] matA matrix A.
	 * @param[in] matB matrix B.
	 * @param[out] matC matrix value: C = A*B.
	 */
	void BLAS_GEMM(CMatrixAS3<as3double> &matA,
			           CMatrixAS3<as3double> &matB,
						     CMatrixAS3<as3double> &matC);

	/*!
	 * @brief Function that computes a matrix-matrix multiplication manually.
	 *
	 * @param[in] matA matrix A.
	 * @param[in] matB matrix B.
	 * @param[out] matC matrix value: C = A*B.
	 */
	void MatrixMatrixMult(CMatrixAS3<as3double> &matA,
			                  CMatrixAS3<as3double> &matB,
						            CMatrixAS3<as3double> &matC);


	/*!
	 * @brief Function that computes the transpose of an AS3-type matrix.
	 *
	 * @param[in] A input matrix A.
	 * @param[out] At transpose of matrix A.
	 */
	void TransposeAS3Matrix(CMatrixAS3<as3double> &A, CMatrixAS3<as3double> &At);

	/*!
	 * @brief Function the creates an identity matrix of dimension n and type T.
	 *
	 * @param[in] n dimension of the (square) matrix.
	 *
	 * @return identity matrix of size n-by-n.
	 */
	template<typename T>
	CMatrixAS3<T> CreateIdentityMatrix(size_t n);


}


// Definitions of the templated functions.
#include "linear_algebra.inl"

