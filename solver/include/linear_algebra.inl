

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
