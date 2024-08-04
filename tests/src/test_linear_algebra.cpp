#include "test_linear_algebra.hpp"


//-----------------------------------------------------------------------------------
// NTest_LA namespace functions.
//-----------------------------------------------------------------------------------


void NTest_LA::InitMatrixSeq
(
 CMatrixAS3<as3double> &A
)
	/*
	 * Generate a sequential set of values in a matrix.
	 */
{
	int s = 0;
	for(size_t i=0; i<A.row(); i++)
		for(size_t j=0; j<A.col(); j++)
			A(i,j) = s++;
}

//-----------------------------------------------------------------------------------

void NTest_LA::InitMatrixRand
(
 CMatrixAS3<as3double> &A
)
	/*
	 * Generate a random set of values in a matrix.
	 */
{
	std::uniform_real_distribution<as3double> dist(-10,50);	
	std::default_random_engine re;
		
	for(size_t i=0; i<A.size(); i++) A[i] = dist(re);
}

//-----------------------------------------------------------------------------------

void NTest_LA::MatrixTranspose
(
 CMatrixAS3<as3double> &A,
 CMatrixAS3<as3double> &B
)
	/*
	 * Compute the correct matrix transpose: B = A^T.
	 */
{
	const size_t m = A.row();
	const size_t n = A.col();
	
	B.resize(n,m);
	for(size_t i=0; i<m; i++)
		for(size_t j=0; j<n; j++)
			B(j,i) = A(i,j);
}

//-----------------------------------------------------------------------------------

void NTest_LA::MatrixMult
(
 CMatrixAS3<as3double> &A,
 CMatrixAS3<as3double> &B,
 CMatrixAS3<as3double> &C
)
	/*
	 * Compute the correct matrix multiplication: C = A*B.
	 */
{
	const size_t m = A.row();
	const size_t n = B.col();
	const size_t c = A.col();

	if( A.col() != B.row() ) FAIL();

	C.resize(m,n);
	for(size_t i=0; i<m; i++)
		for(size_t j=0; j<n; j++)
			for(size_t k=0; k<c; k++)
				C(i,j) += A(i,k)*B(k,j);
}

//-----------------------------------------------------------------------------------

void NTest_LA::MatrixErrorNormLinf
(
 CMatrixAS3<as3double> &A,
 CMatrixAS3<as3double> &B,
 as3double              tol
)
	/*
	 * Check the L1 norm of the difference in matrix: A and B.
	 */
{
	// Numeric minimum tolerance. The 10 is a safety factor.
	const as3double mintol = 10.0*std::numeric_limits<as3double>::epsilon();

	if( A.row() != B.row() ) FAIL();
	if( A.col() != B.col() ) FAIL();

	as3double errmax = 0.0;
	for(size_t i=0; i<A.size(); i++) 
	{		
		const as3double reltol = mintol*std::max( std::fabs(A[i]), std::fabs(B[i]) );
		const as3double   diff = std::fabs(A[i] - B[i]);

		if( diff > reltol ) errmax = std::max( diff, errmax );
	}

	ASSERT_TRUE( errmax < tol );
}

//-----------------------------------------------------------------------------------

void NTest_LA::InverseMatrix3x3
(
 CMatrixAS3<as3double> &A,
 CMatrixAS3<as3double> &B
)
	/*
	 * Computes the analytical inverse of a 3x3 matrix: B = inv(A).
	 */
{
	as3double zero = 0.0;
	
	A.resize(3,3);
	B.resize(3,3);

	// Generate a matrix with a 'nice' condition number.
	A(0,0) = 2.0; A(0,1) = 1.0; A(0,2) = 3.0;
	A(1,0) = 0.0; A(1,1) = 2.0; A(1,2) = 4.0;
	A(2,0) = 1.0; A(2,1) = 1.0; A(2,2) = 2.0;


	as3double det = +A(0,0)*( A(1,1)*A(2,2) - A(2,1)*A(1,2) )
	                -A(0,1)*( A(1,0)*A(2,2) - A(1,2)*A(2,0) )
	                +A(0,2)*( A(1,0)*A(2,1) - A(1,1)*A(2,0) );

	B(0,0) =  ( A(1,1)*A(2,2) - A(2,1)*A(1,2) )/det;
	B(0,1) = -( A(0,1)*A(2,2) - A(0,2)*A(2,1) )/det;
	B(0,2) =  ( A(0,1)*A(1,2) - A(0,2)*A(1,1) )/det;
	B(1,0) = -( A(1,0)*A(2,2) - A(1,2)*A(2,0) )/det;
	B(1,1) =  ( A(0,0)*A(2,2) - A(0,2)*A(2,0) )/det;
	B(1,2) = -( A(0,0)*A(1,2) - A(1,0)*A(0,2) )/det;
	B(2,0) =  ( A(1,0)*A(2,1) - A(2,0)*A(1,1) )/det;
	B(2,1) = -( A(0,0)*A(2,1) - A(2,0)*A(0,1) )/det;
	B(2,2) =  ( A(0,0)*A(1,1) - A(1,0)*A(0,1) )/det;
}

//-----------------------------------------------------------------------------------

void NTest_LA::PrintAS3Matrix
(
 CMatrixAS3<as3double> &A,
 std::string            g,
 unsigned short         p
)
	/*
	 * Prints an AS3-type matrix, for debugging.
	 */
{	
	std::cout << "Matrix: " << g << std::endl;		
	for(size_t i=0; i<A.row(); i++)		
	{
		for(size_t j=0; j<A.col(); j++)
			std::cout << std::showpos << std::scientific << std::setprecision(p) << std::setw(p) << A(i,j) << ", ";
		std::cout << std::endl;
	}
}
