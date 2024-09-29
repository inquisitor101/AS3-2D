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
	 * @param[in,out] A inverse of the matrix A (over-written). 
	 */
	void InverseMatrix(CMatrixAS3<as3double> &A);

	/*!
	 * @brief Function that computes a generalized matrix-matrix multiplication via BLAS (double-precision).
	 *
	 * @param[in] matA matrix A.
	 * @param[in] matB matrix B.
	 *
	 * return matrix value: C = A*B.
	 */
	CMatrixAS3<as3double> BLAS_GEMM(CMatrixAS3<as3double> &matA,
			                            CMatrixAS3<as3double> &matB);

	/*!
	 * @brief Function that computes a matrix-matrix multiplication manually.
	 *
	 * @param[in] matA matrix A.
	 * @param[in] matB matrix B.
	 *
	 * @return matrix value: C = A*B.
	 */
	CMatrixAS3<as3double> MatrixMatrixMult(CMatrixAS3<as3double> &matA,
			                                   CMatrixAS3<as3double> &matB);

	/*!
	 * @brief Function that computes a matrix-vector multiplication manually.
	 *
	 * @param[in] matA matrix A.
	 * @param[in] vecB vector B.
	 *
	 * @return vector value: {C} = A*{B}, where {} can be multiple vectors.
	 */
	CMatrixAS3<as3double> MatrixVectorMult(CMatrixAS3<as3double> &matA,
			                                   CMatrixAS3<as3double> &vecB);

	/*!
	 * @brief Function that computes a matrix-vector multiplication manually, with vecB being transposed.
	 *
	 * @param[in] matA matrix A.
	 * @param[in] vecB vector B.
	 * @param[out] vecC vector: {C} = A*{B}, where {} can be multiple vectors.
	 */
	void MatrixVectorTransMult(CMatrixAS3<as3double> &matA,
			                       CMatrixAS3<as3double> &vecB,
												     CMatrixAS3<as3double> &vecC);

	/*!
	 * @brief Function that computes the transpose of an AS3-type matrix.
	 *
	 * @param[in] A input matrix A.
	 * @param[out] At transpose of matrix A.
	 */
	void TransposeAS3Matrix(CMatrixAS3<as3double> &A, CMatrixAS3<as3double> &At);

	/*!
	 * @brief Function that creates the output of the transpose of an AS3-type matrix.
	 *
	 * @param[in] A input matrix A.
	 *
	 * @return the transpose of matrix A.
	 */
	template<typename T>
	CMatrixAS3<T> TransposeAS3Matrix(CMatrixAS3<T> &A);

	/*!
	 * @brief Function the creates an identity matrix of dimension n and type T.
	 *
	 * @param[in] n dimension of the (square) matrix.
	 *
	 * @return identity matrix of size n-by-n.
	 */
	template<typename T>
	CMatrixAS3<T> CreateIdentityMatrix(size_t n);

	/*!
	 * @brief Function the creates the output of a kronecker (outer-)product of two matrices.
	 *
	 * @param[in] A reference to the left-hand matrix.
	 * @param[in] B reference to the right-hand matrix.
	 *
	 * @return outer-product matrix: C = kron(A,B).
	 */
	template<typename T>
	CMatrixAS3<T> KroneckerProduct(CMatrixAS3<T> &A, 
			                           CMatrixAS3<T> &B);
}


// Definitions of the templated functions.
#include "linear_algebra.inl"

