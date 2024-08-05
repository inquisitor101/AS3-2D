#include "linear_algebra.hpp"



/*!
 * @brief Function prototypes for the Lapack routines used.
 */
extern "C"
{
	// LU decomoposition of a general matrix (double-precision).
  void dgetrf_(int*, int *, double*, int*, int*, int*);

  // Generate inverse of a matrix given its LU decomposition (double-precision).
  void dgetri_(int*, double*, int*, int*, double*, int*, int*);

	// Compute a generalized matrix-matrix multiplication (double-precision).
	int dgemm_(char*, char*, int*, int*, int*, double*, double*, 
			       int*, double*, int*, double*, double*, int*);
}



//-----------------------------------------------------------------------------------
// NLinearAlgebra namespace functions.
//-----------------------------------------------------------------------------------


void NLinearAlgebra::InverseMatrix
(
 CMatrixAS3<as3double> &A
)
 /*
	* Function that computes the inverse of the given matrix (double-precision).
	*/
{
	// Check the matrix is indeed a square matrix.
	if( A.row() != A.col() ) ERROR("Cannot invert a non-square matrix.");

	// Cast to integer safely.
	int n = static_cast<int>( A.row() );
  
	// Allocate the memory for the work arrays.
  as3vector1d<int>       ipivVec(n);
  as3vector1d<as3double> workVec(n);

  // Determine the pointers for the data in the vectors.
  int *ipiv = ipivVec.data();

  as3double *work = workVec.data();
  as3double *AA   = A.data();

  // Call the appropriate Lapack functions to compute the inverse.
  int ierr;

	// Note, this assumes double-precision.
  dgetrf_(&n, &n, AA, &n, ipiv, &ierr);
  if(ierr != 0) ERROR("Matrix is singular");

  dgetri_(&n, AA, &n, ipiv, work, &n, &ierr);
  if(ierr != 0) ERROR("Matrix inversion failed");
}

//-----------------------------------------------------------------------------------

void NLinearAlgebra::BLAS_GEMM
(
 CMatrixAS3<as3double> &matA,
 CMatrixAS3<as3double> &matB,
 CMatrixAS3<as3double> &matC
)
 /*
	* Function that computes a generalized matrix-matrix multiplication via BLAS (double-precision).
	* Note, this routine allocates the memory for matC.
	*/
{
	// Note, the GEMM BLAS routine assumes column-major, since the implementation is in Fortran.
	// The CMatrixAS3 uses row-major, in line with the C++ data structure format. To avoid this 
	// issue, notice no transpose is computed in the DGEMM function, but the matrices have been
	// switched, given the fact: C^T = B^T * A^T. 
	// Here, the returned matrix C^T is in column-major, which is equivalent to C, since we use
	// row-major. For simplicity, the dimensions of the matrices are:
	// 
	// A: Ma*Na, B: Mb*Nb  <==> A^T: Na*Ma, B^T: Nb*Mb.
	// 
	// Objective: C = A*B.
	//
	// We achieve this as follows,
	// C^T = B^T * A^T, where C^T: Nb*Ma.
	// ... in the following, 
	//  m = Nb,
	//  n = Ma,
	//  c = Mb (common dimension)
	//---------------------------

	// Ensure the common dimension of the matrices is correct.
	if( matA.col() != matB.row() ) ERROR("Matrix dimensions are incorrect.");

	char tA = 'N';
	char tB = 'N';

	double one  = 1.0;
	double zero = 0.0;

	int m = static_cast<int>(matB.col()); // rows of  left-matrix, i.e.: B^T.
	int n = static_cast<int>(matA.row()); // cols of right-matrix, i.e.: A^T.
	int c = static_cast<int>(matB.row()); // common size.

	// Allocate the matrix to the correct size in row-major form, i.e. rows(A)*col(B).
	matC.resize(n, m);

	as3double *A = matB.data(); // switched to B.
	as3double *B = matA.data(); // switched to A.
	as3double *C = matC.data();

	// Call the generalized matrix-matrix routine (double-precision).
	dgemm_(&tA, &tB, &m, &n, &c, &one, A, &m, B, &c, &zero, C, &m);
}

//-----------------------------------------------------------------------------------

void NLinearAlgebra::MatrixMatrixMult
(
 CMatrixAS3<as3double> &matA,
 CMatrixAS3<as3double> &matB,
 CMatrixAS3<as3double> &matC
)
 /*
	* Function that computes a matrix-matrix multiplication manually.
	* Note, this routine allocates the memory for matC.
	*/
{
	// Ensure the common dimension of the matrices is correct.
	if( matA.col() != matB.row() ) ERROR("Matrix dimensions are incorrect.");

	// Extract the numbner of rows and columns of C and the common dimension.
	const size_t nRow = matA.row();
	const size_t nCol = matB.col();
	const size_t nVec = matA.col(); // Common dimension.

	// Allocate the matrix to the correct size in row-major form, i.e. rows(A)*col(B).
	matC.resize(nRow, nCol);

	// Perform actual matrix-matrix multiplication.
	for(size_t i=0; i<nRow; i++)
	{
		for(size_t j=0; j<nCol; j++)
		{
			matC(i,j) = 0.0; 
			for(size_t k=0; k<nVec; k++)
			{
				matC(i,j) += matA(i,k)*matB(k,j);
			}
		}
	}
}

//-----------------------------------------------------------------------------------

void NLinearAlgebra::TransposeAS3Matrix
(
 CMatrixAS3<as3double> &A,
 CMatrixAS3<as3double> &At
)
 /*
	* Function that computes the transpose of an AS3-type matrix.
	* Note, this routine allocates the memory for At.
	*/
{
	// Extract the dimensions.
	const size_t m = A.row();
	const size_t n = A.col();

	// Allocate memory for the transposed matrix.
	At.resize( n, m );

	// Transpose the matrix.
	for(size_t i=0; i<m; i++)
	{
		for(size_t j=0; j<n; j++)
		{
			At(j,i) = A(i,j);
		}
	}
}





