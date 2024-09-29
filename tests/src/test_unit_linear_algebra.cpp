#include "test_unit_linear_algebra.hpp"




namespace 
{

	TEST(LinearAlgebra, Transpose) 
	{
		const size_t m = 10;
		const size_t n = 15;
		CMatrixAS3<as3double> A(m,n);
		CMatrixAS3<as3double> B;
		CMatrixAS3<as3double> C;

		NTest_LA::InitMatrixSeq(A);
		NTest_LA::MatrixTranspose(A,B);
		
		NLinearAlgebra::TransposeAS3Matrix(A, C);

		if( B.row() != C.row() ) FAIL();
		if( B.col() != C.col() ) FAIL();

		NTest_LA::MatrixErrorNormLinf(B, C, 1.0e-14);
	}


	TEST(LinearAlgebra, MatrixMult) 
	{
		const size_t m = 10;
		const size_t n = 15;
		const size_t c = 20;
		
		CMatrixAS3<as3double> A(m,c);
		CMatrixAS3<as3double> B(c,n);	

		NTest_LA::InitMatrixRand(A);
		NTest_LA::InitMatrixRand(B);

		CMatrixAS3<as3double> E;

		CMatrixAS3<as3double> C = NLinearAlgebra::MatrixMatrixMult(A, B);
		CMatrixAS3<as3double> D = NLinearAlgebra::BLAS_GEMM(A, B);
		NTest_LA::MatrixMult(A, B, E);
		
		NTest_LA::MatrixErrorNormLinf(C, E, 1.0e-14);
		NTest_LA::MatrixErrorNormLinf(D, E, 1.0e-14);
	}


	TEST(LinearAlgebra, MatrixInv)
	{
		CMatrixAS3<as3double> A;
		CMatrixAS3<as3double> B;

		NTest_LA::InverseMatrix3x3(A, B);
		
		CMatrixAS3<as3double> C = A;
		NLinearAlgebra::InverseMatrix(C);
		
		NTest_LA::MatrixErrorNormLinf(B, C, 1.0e-14);
	
		// Create the identity matrix.
		CMatrixAS3<as3double> I(10,10);
		for(size_t k=0; k<I.row(); k++) I(k,k) = 1.0;

		A.resize(I.row(), I.col());
		for(size_t k=0; k<A.row(); k++) A(k,k) = (as3double) k+1.0;

		C = A;
		NLinearAlgebra::InverseMatrix(C);

		CMatrixAS3<as3double> D = NLinearAlgebra::BLAS_GEMM(C, A);
		NTest_LA::MatrixErrorNormLinf(D, I, 1.0e-14);

		CMatrixAS3<as3double> E = NLinearAlgebra::MatrixMatrixMult(C, A);
		NTest_LA::MatrixErrorNormLinf(E, I, 1.0e-14);
	}
}





