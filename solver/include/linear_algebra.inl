/*!
 * @file
 * @brief Inline implementations for dense matrix operations and tensor-product matrix utilities.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the templated function in NLinearAlgebra.
//-----------------------------------------------------------------------------------

/*!
 * @brief Function that creates an identity matrix of dimension n and type T.
 *
 * @tparam T matrix scalar type.
 *
 * @param[in] n dimension of the square identity matrix.
 *
 * @return identity matrix of size n-by-n.
 */
template<typename T>
CMatrixAS3<T> NLinearAlgebra::CreateIdentityMatrix
(
 std::size_t n
)
 /*
	* Function that creates an identity matrix and returns it via RVO.
	*/
{
	// Allocate the matrix.
	CMatrixAS3<T> matrix(n,n);

	// Avoid floating-precision, by casting the data types once.
	const T one = static_cast<T>( 1.0 );

	// Specify ones on its diagonal.
	for(std::size_t i=0; i<n; i++) matrix(i,i) = one;

	// Return the matrix via copy elision (RVO).
	return matrix;
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that computes the Kronecker product of two matrices.
 *
 * @tparam T matrix scalar type.
 *
 * @param[in] A input matrix.
 * @param[in] B right input matrix.
 *
 * @return Kronecker product of A and B.
 */
template<typename T>
CMatrixAS3<T> NLinearAlgebra::KroneckerProduct
(
 CMatrixAS3<T> &A,
 CMatrixAS3<T> &B
)
 /*
	* Function that computes the kronecker product of matrix A and B and returns it via RVO.
	*/
{
	// Obtain the dimensions of the input matrices 
	const std::size_t ma = A.row();
	const std::size_t na = A.col();
	const std::size_t mb = B.row();
	const std::size_t nb = B.col();

	// Deduce the dimensions of the outer-product matrix C = kron(A,B).
	const std::size_t mc = ma*mb;
	const std::size_t nc = na*nb;

	// Allocate the matrix.
	CMatrixAS3<T> C(mc,nc);

	// Perform the actual outer-product.
	for(std::size_t ia=0; ia<ma; ia++)
	{
		for(std::size_t ib=0; ib<mb; ib++)
		{
			for(std::size_t ja=0; ja<na; ja++)
			{
				for(std::size_t jb=0; jb<nb; jb++)
				{
					C(ia*mb+ib, ja*nb+jb) = A(ia,ja)*B(ib,jb);
				}
			}
		}
	}

	// Return the matrix C via copy elision (RVO).
	return C;
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that creates the output of the transpose of an AS3-type matrix.
 *
 * @tparam T matrix scalar type.
 *
 * @param[in] A input matrix.
 *
 * @return transpose of A.
 */
template<typename T>
CMatrixAS3<T> NLinearAlgebra::TransposeAS3Matrix
(
 CMatrixAS3<T> &A
)
 /*
	* Function that computes the transpose of a generic AS3 matrix and returns it via RVO.
	*/
{
	// Extract the dimensions.
	const std::size_t m = A.row();
	const std::size_t n = A.col();

	// Allocate memory for the transposed matrix.
	CMatrixAS3<T> At( n, m );

	// Transpose the matrix.
	for(std::size_t i=0; i<m; i++)
	{
		for(std::size_t j=0; j<n; j++)
		{
			At(j,i) = A(i,j);
		}
	}

	// Return the matrix At via copy elision (RVO).
	return At;
}





