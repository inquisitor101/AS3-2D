

//-----------------------------------------------------------------------------------
// Implementation of the templated function in NLinearAlgebra.
//-----------------------------------------------------------------------------------

template<typename T>
CMatrixAS3<T> NLinearAlgebra::CreateIdentityMatrix
(
 size_t n
)
 /*
	* Function that creates an identity matrix and returns it via RVO.
	*/
{
	// Allocate the matrix.
	CMatrixAS3<T> matrix(n,n);

	// Avoid floating-precision, by casting the data types onces.
	const T one = static_cast<T>( 1.0 );

	// Specify ones on its diagonal.
	for(size_t i=0; i<n; i++) matrix(i,i) = one;

	// Return the matrix via copy elision (RVO).
	return matrix;
}

//-----------------------------------------------------------------------------------

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
	const size_t ma = A.row();
	const size_t na = A.col();
	const size_t mb = B.row();
	const size_t nb = B.col();

	// Deduce the dimensions of the outer-product matrix C = kron(A,B).
	const size_t mc = ma*mb;
	const size_t nc = na*nb;

	// Allocate the matrix.
	CMatrixAS3<T> C(mc,nc);

	// Perform the actual outer-product.
	for(size_t ia=0; ia<ma; ia++)
	{
		for(size_t ib=0; ib<mb; ib++)
		{
			for(size_t ja=0; ja<na; ja++)
			{
				for(size_t jb=0; jb<nb; jb++)
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
	const size_t m = A.row();
	const size_t n = A.col();

	// Allocate memory for the transposed matrix.
	CMatrixAS3<T> At( n, m );

	// Transpose the matrix.
	for(size_t i=0; i<m; i++)
	{
		for(size_t j=0; j<n; j++)
		{
			At(j,i) = A(i,j);
		}
	}

	// Return the matrix At via copy elision (RVO).
	return At;
}





