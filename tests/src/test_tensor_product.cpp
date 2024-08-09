#include "test_tensor_product.hpp"


//-----------------------------------------------------------------------------------
// CTest_TP member functions.
//-----------------------------------------------------------------------------------


void CTest_TP::SetUp
(
 void
)
 /*
	* Function that serves as a constructor for the object.
	*/
{
	// Suppress cout
	auto buffer = std::cout.rdbuf(nullptr);

	mConfig  = std::make_unique<CConfig>("param_t1.cfg");
	ASSERT_GT(mConfig->GetnZone(), 0);

	mElement.resize(mConfig->GetnZone());
	for(unsigned short iZone=0; iZone<mConfig->GetnZone(); iZone++) 
		mElement[iZone] = std::make_unique<CStandardElement>(mConfig.get(), iZone);

	// Restore cout.
	std::cout.rdbuf(buffer);
}

//-----------------------------------------------------------------------------------

void CTest_TP::KroneckerProduct
(
 CMatrixAS3<as3double> &A,
 CMatrixAS3<as3double> &B,
 CMatrixAS3<as3double> &C
)
 /*
	* Function that computes the kronecker/outer product of 2 matrices, as such: C = kron(A,B).
	*/
{
	const size_t ma = A.row();
	const size_t na = A.col();
	const size_t mb = B.row();
	const size_t nb = B.col();

	const size_t mc = ma*mb;
	const size_t nc = na*nb;

	C.resize(mc, nc);

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
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckVolumeSol
(
 CStandardElement *element
)
 /*
	* Function that checks the volume tensor-product interpolation.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;
	
	CMatrixAS3<as3double> A;   // input solution.
	CMatrixAS3<as3double> At;  // its tranpose.
	CMatrixAS3<as3double> B;   // benchmark result.
	CMatrixAS3<as3double> Bt;  // benchmark result transposed.
	CMatrixAS3<as3double> C;   // tensor-product result.
	CMatrixAS3<as3double> L2D; // interpolation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> L1D = element->GetLagrangeInt1D();
	KroneckerProduct(L1D, L1D, L2D);
	NTest_LA::MatrixMult(L2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m2);
	
	// Compile-time tensor-product implementation.
	tensor->Volume(c, A.data(), C.data(), nullptr, nullptr);
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomVolume(n, c, m,
			                 element->GetLagrangeInt1DTrans().data(), 
											 nullptr, 
											 A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckVolumeDerSolR
(
 CStandardElement *element
)
 /*
	* Function that checks the volume tensor-product r-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;
	
	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> drL2D; // r-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	KroneckerProduct(L1D, dL1D, drL2D);
	NTest_LA::MatrixMult(drL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m2);

	// Compile-time tensor product implementation.
	tensor->Volume(c, A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomVolume(n, c, m,
			                 element->GetLagrangeInt1DTrans().data(), 
											 element->GetDerLagrangeInt1DTrans().data(), 
											 A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckVolumeDerSolS
(
 CStandardElement *element
)
 /*
	* Function that checks the volume tensor-product s-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> dsL2D; // s-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	KroneckerProduct(dL1D, L1D, dsL2D);
	NTest_LA::MatrixMult(dsL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m2);

	// Compile-time tensor product implementation.
	tensor->Volume(c, A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomVolume(n, c, m,
			                 element->GetLagrangeInt1DTrans().data(), 
											 element->GetDerLagrangeInt1DTrans().data(), 
											 A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurfaceSolIMIN
(
 CStandardElement *element
)
 /*
	* Function that checks the IMIN surface tensor-product interpolation.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;   // input solution.
	CMatrixAS3<as3double> At;  // its transpose.
	CMatrixAS3<as3double> B;   // benchmark result.
	CMatrixAS3<as3double> Bt;  // benchmark result transposed.
	CMatrixAS3<as3double> C;   // tensor-product result.
	CMatrixAS3<as3double> L2D; // interpolation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> L1D = element->GetLagrangeInt1D();
	
	CMatrixAS3<as3double> Limin(1, n);
	Limin[0] = 1.0;
	
	KroneckerProduct(L1D, Limin, L2D);
	NTest_LA::MatrixMult(L2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceIMIN(c, A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceIMIN(n, c, m,
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMinFace1D().data(),
											      A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurfaceSolIMAX
(
 CStandardElement *element
)
 /*
	* Function that checks the IMAX surface tensor-product interpolation.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;   // input solution.
	CMatrixAS3<as3double> At;  // its transpose.
	CMatrixAS3<as3double> B;   // benchmark result.
	CMatrixAS3<as3double> Bt;  // benchmark result transposed.
	CMatrixAS3<as3double> C;   // tensor-product result.
	CMatrixAS3<as3double> L2D; // interpolation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> L1D = element->GetLagrangeInt1D();
	
	CMatrixAS3<as3double> Limax(1, n);
	Limax[n-1] = 1.0;
	
	KroneckerProduct(L1D, Limax, L2D);
	NTest_LA::MatrixMult(L2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceIMAX(c, A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceIMAX(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      nullptr,
														nullptr,
											      A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurfaceSolJMIN
(
 CStandardElement *element
)
 /*
	* Function that checks the JMIN surface tensor-product interpolation.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;   // input solution.
	CMatrixAS3<as3double> At;  // its transpose.
	CMatrixAS3<as3double> B;   // benchmark result.
	CMatrixAS3<as3double> Bt;  // benchmark result transposed.
	CMatrixAS3<as3double> C;   // tensor-product result.
	CMatrixAS3<as3double> L2D; // interpolation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> L1D = element->GetLagrangeInt1D();

	// Compute the surface values.
	L2D.resize(m, n2);
	for(size_t i=0; i<L1D.row(); i++)
		for(size_t j=0; j<L1D.col(); j++)
			L2D(i,j) = L1D(i,j);

	NTest_LA::MatrixMult(L2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceJMIN(c, A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceJMIN(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      nullptr,
														nullptr,
											      A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurfaceSolJMAX
(
 CStandardElement *element
)
 /*
	* Function that checks the JMAX surface tensor-product interpolation.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;   // input solution.
	CMatrixAS3<as3double> At;  // its transpose.
	CMatrixAS3<as3double> B;   // benchmark result.
	CMatrixAS3<as3double> Bt;  // benchmark result transposed.
	CMatrixAS3<as3double> C;   // tensor-product result.
	CMatrixAS3<as3double> L2D; // interpolation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> L1D = element->GetLagrangeInt1D();

	// Compute the surface values.
	L2D.resize(m, n2);
	size_t s = n2 - n;
	for(size_t i=0; i<L1D.row(); i++)
		for(size_t j=0; j<L1D.col(); j++)
			L2D(i,j+s) = L1D(i,j);

	NTest_LA::MatrixMult(L2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceJMAX(c, A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceJMAX(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      nullptr,
														nullptr,
											      A.data(), C.data(), nullptr, nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDrIMIN
(
 CStandardElement *element
)
 /*
	* Function that checks the IMIN surface tensor-product r-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> drL2D; // r-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeSol1D();
	
	CMatrixAS3<as3double> dLimin(1, n);	
	for(size_t i=0; i<dL1D.col(); i++) dLimin[i] = dL1D(0, i);

	KroneckerProduct(L1D, dLimin, drL2D);
	NTest_LA::MatrixMult(drL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceIMIN(c, A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceIMIN(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMinFace1D().data(),
											      A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDrIMAX
(
 CStandardElement *element
)
 /*
	* Function that checks the IMAX surface tensor-product r-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> drL2D; // r-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeSol1D();
	
	CMatrixAS3<as3double> dLimax(1, n);	
	for(size_t i=0; i<n; i++) dLimax[i] = dL1D(n-1, i);

	KroneckerProduct(L1D, dLimax, drL2D);
	NTest_LA::MatrixMult(drL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceIMAX(c, A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceIMAX(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMaxFace1D().data(),
											      A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDrJMIN
(
 CStandardElement *element
)
 /*
	* Function that checks the JMIN surface tensor-product r-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> drL2D; // r-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	
	// Compute the surface values.
	drL2D.resize(m, n2);
	for(size_t i=0; i<dL1D.row(); i++)
		for(size_t j=0; j<dL1D.col(); j++)
			drL2D(i,j) = dL1D(i,j);

	//KroneckerProduct(L1D, dLjmin, drL2D);
	NTest_LA::MatrixMult(drL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceJMIN(c, A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceJMIN(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMinFace1D().data(),
											      A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDrJMAX
(
 CStandardElement *element
)
 /*
	* Function that checks the JMAX surface tensor-product r-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> drL2D; // r-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	
	// Compute the surface values.
	drL2D.resize(m, n2);
	size_t s = n2 - n;
	for(size_t i=0; i<dL1D.row(); i++)
		for(size_t j=0; j<dL1D.col(); j++)
			drL2D(i,j+s) = dL1D(i,j);

	//KroneckerProduct(L1D, dLjmin, drL2D);
	NTest_LA::MatrixMult(drL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceJMAX(c, A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceJMAX(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMaxFace1D().data(),
											      A.data(), nullptr, C.data(), nullptr);

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDsIMIN
(
 CStandardElement *element
)
 /*
	* Function that checks the IMIN surface tensor-product s-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> dsL2D; // s-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	
	CMatrixAS3<as3double> dLimin(1, n);
	dLimin[0] = 1.0;

	KroneckerProduct(dL1D, dLimin, dsL2D);
	NTest_LA::MatrixMult(dsL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceIMIN(c, A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceIMIN(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMinFace1D().data(),
											      A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDsIMAX
(
 CStandardElement *element
)
 /*
	* Function that checks the IMAX surface tensor-product s-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> dsL2D; // s-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	
	CMatrixAS3<as3double> dLimax(1, n);
	dLimax[n-1] = 1.0;

	KroneckerProduct(dL1D, dLimax, dsL2D);
	NTest_LA::MatrixMult(dsL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceIMAX(c, A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceIMAX(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMaxFace1D().data(),
											      A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDsJMIN
(
 CStandardElement *element
)
 /*
	* Function that checks the JMIN surface tensor-product s-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> dsL2D; // s-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeSol1D();
	
	CMatrixAS3<as3double> dLjmin(1, n);	
	for(size_t i=0; i<n; i++) dLjmin[i] = dL1D(0, i);

	KroneckerProduct(dLjmin, L1D, dsL2D);
	NTest_LA::MatrixMult(dsL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceJMIN(c, A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceJMIN(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMinFace1D().data(),
											      A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckSurface_dSolDsJMAX
(
 CStandardElement *element
)
 /*
	* Function that checks the JMAX surface tensor-product s-derivative.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;

	CMatrixAS3<as3double> A;     // input solution.
	CMatrixAS3<as3double> At;    // its transpose.
	CMatrixAS3<as3double> B;     // benchmark result.
	CMatrixAS3<as3double> Bt;    // benchmark result transposed.
	CMatrixAS3<as3double> C;     // tensor-product result.
	CMatrixAS3<as3double> dsL2D; // s-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, n2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeSol1D();
	
	CMatrixAS3<as3double> dLjmax(1, n);	
	for(size_t i=0; i<n; i++) dLjmax[i] = dL1D(n-1, i);

	KroneckerProduct(dLjmax, L1D, dsL2D);
	NTest_LA::MatrixMult(dsL2D, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, m);

	// Compile-time tensor product implementation.
	tensor->SurfaceJMAX(c, A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);

	// Reset solution.
	C.reset();

	// Run-time tensor product implementation.
	tensor->CustomSurfaceJMAX(n, c, m, 
			                      element->GetLagrangeInt1DTrans().data(), 
											      element->GetDerLagrangeInt1DTrans().data(),
														element->GetDerLagrangeMaxFace1D().data(),
											      A.data(), nullptr, nullptr, C.data());

	// Check error.
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckResidualVolumeSource
(
 CStandardElement *element
)
 /*
	* Function that checks the residual volume tensor-product for source terms.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;
	
	CMatrixAS3<as3double> A;    // input solution.
	CMatrixAS3<as3double> At;   // its tranpose.
	CMatrixAS3<as3double> B;    // benchmark result.
	CMatrixAS3<as3double> Bt;   // benchmark result transposed.
	CMatrixAS3<as3double> C;    // tensor-product result.
	CMatrixAS3<as3double> L2Dt; // transpose of interpolation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, m2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double> L1Dt = element->GetLagrangeInt1DTrans();
	KroneckerProduct(L1Dt, L1Dt, L2Dt);
	NTest_LA::MatrixMult(L2Dt, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, n2);
	
	// Compile-time tensor-product implementation.
	tensor->ResidualVolume(c, A.data(), nullptr, nullptr, C.data());
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckResidualVolumeDerSolR
(
 CStandardElement *element
)
 /*
	* Function that checks the residual volume tensor-product for r-differentiated terms.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;
	
	CMatrixAS3<as3double> A;      // input solution.
	CMatrixAS3<as3double> At;     // its tranpose.
	CMatrixAS3<as3double> B;      // benchmark result.
	CMatrixAS3<as3double> Bt;     // benchmark result transposed.
	CMatrixAS3<as3double> C;      // tensor-product result.
	CMatrixAS3<as3double> drL2D;  // r-differentiation matrix in 2D.
	CMatrixAS3<as3double> drL2Dt; // transpose of r-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, m2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	KroneckerProduct(L1D, dL1D, drL2D);
	NTest_LA::MatrixTranspose(drL2D, drL2Dt);

	NTest_LA::MatrixMult(drL2Dt, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, n2);
	
	// Compile-time tensor-product implementation.
	tensor->ResidualVolume(c, nullptr, A.data(), nullptr, C.data());
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckResidualVolumeDerSolS
(
 CStandardElement *element
)
 /*
	* Function that checks the residual volume tensor-product for s-differentiated terms.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;
	
	CMatrixAS3<as3double> A;      // input solution.
	CMatrixAS3<as3double> At;     // its tranpose.
	CMatrixAS3<as3double> B;      // benchmark result.
	CMatrixAS3<as3double> Bt;     // benchmark result transposed.
	CMatrixAS3<as3double> C;      // tensor-product result.
	CMatrixAS3<as3double> dsL2D;  // s-differentiation matrix in 2D.
	CMatrixAS3<as3double> dsL2Dt; // transpose of s-differentiation matrix in 2D.

	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A.resize(c, m2);
	NTest_LA::InitMatrixRand(A);
	NTest_LA::MatrixTranspose(A, At);

	// Form the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();
	KroneckerProduct(dL1D, L1D, dsL2D);
	NTest_LA::MatrixTranspose(dsL2D, dsL2Dt);

	NTest_LA::MatrixMult(dsL2Dt, At, B);
	NTest_LA::MatrixTranspose(B, Bt);

	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, n2);
	
	// Compile-time tensor-product implementation.
	tensor->ResidualVolume(c, nullptr, nullptr, A.data(), C.data());
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}

//-----------------------------------------------------------------------------------

void CTest_TP::CheckResidualVolumeTotal
(
 CStandardElement *element
)
 /*
	* Function that checks the residual volume tensor-product for source and rs-differentiated terms.
	*/
{
	// Relative error tolerance.
	const as3double tol = 1.0e-12;
	
	CMatrixAS3<as3double> A1;     // input solution: source terms.
	CMatrixAS3<as3double> A2;     // input solution: r-diff terms.
	CMatrixAS3<as3double> A3;     // input solution: s-diff terms.
	CMatrixAS3<as3double> At1;    // input solution tranpose: source terms.
	CMatrixAS3<as3double> At2;    // input solution tranpose: r-diff terms.
	CMatrixAS3<as3double> At3;    // input solution tranpose: s-diff terms.
	CMatrixAS3<as3double> B;      // benchmark result.
	CMatrixAS3<as3double> Bt;     // benchmark result transposed: containing all terms.
	CMatrixAS3<as3double> Bt1;    // benchmark result transposed: source terms.
	CMatrixAS3<as3double> Bt2;    // benchmark result transposed: r-diff terms.
	CMatrixAS3<as3double> Bt3;    // benchmark result transposed: s-diff terms.
	CMatrixAS3<as3double> C;      // tensor-product result.
	CMatrixAS3<as3double> L2Dt;   // transpose of interpolation matrix in 2D.
	CMatrixAS3<as3double> drL2D;  // r-differentiation matrix in 2D.
	CMatrixAS3<as3double> drL2Dt; // transpose of r-differentiation matrix in 2D.
	CMatrixAS3<as3double> dsL2D;  // s-differentiation matrix in 2D.
	CMatrixAS3<as3double> dsL2Dt; // transpose of s-differentiation matrix in 2D.


	// Number of integration and solution points in 1D.
	const size_t c  = 3;
	const size_t n  = element->GetnSol1D();
	const size_t m  = element->GetnInt1D();
	const size_t n2 = element->GetnSol2D();
	const size_t m2 = element->GetnInt2D();

	// Initialize the matrix A randomly.
	A1.resize(c, m2);
	NTest_LA::InitMatrixRand(A1);
	NTest_LA::MatrixTranspose(A1, At1);

	// Initialize arbitrarily its transpose for the rs-differentiation terms.
	A2 = A1; NTest_LA::InitMatrixRand(A2, -5.0, 5.0);
	A3 = A1; NTest_LA::InitMatrixRand(A3,  0.0, 3.0);
	NTest_LA::MatrixTranspose(A2, At2);
	NTest_LA::MatrixTranspose(A3, At3);


	// Required basis functions for the benchmark values.
	CMatrixAS3<as3double>  L1D = element->GetLagrangeInt1D();
	CMatrixAS3<as3double> L1Dt = element->GetLagrangeInt1DTrans();
	CMatrixAS3<as3double> dL1D = element->GetDerLagrangeInt1D();

	// Form the benchmark values for the source terms.
	KroneckerProduct(L1Dt, L1Dt, L2Dt);
	NTest_LA::MatrixMult(L2Dt, At1, B);
	NTest_LA::MatrixTranspose(B, Bt1);

	// Form the benchmark values for the r-differentiated terms.
	KroneckerProduct(L1D, dL1D, drL2D);
	NTest_LA::MatrixTranspose(drL2D, drL2Dt);
	NTest_LA::MatrixMult(drL2Dt, At2, B);
	NTest_LA::MatrixTranspose(B, Bt2);

	// Form the benchmark values for the s-differentiated terms.
	KroneckerProduct(dL1D, L1D, dsL2D);
	NTest_LA::MatrixTranspose(dsL2D, dsL2Dt);
	NTest_LA::MatrixMult(dsL2Dt, At3, B);
	NTest_LA::MatrixTranspose(B, Bt3);

	// Accumulate all the term contributions in Bt.
	Bt.resize( Bt1.row(), Bt1.col() );
	for(size_t i=0; i<Bt.size(); i++) Bt[i] = Bt1[i] + Bt2[i] + Bt3[i];


	// Create a tensor product object.
	std::unique_ptr<ITensorProduct> tensor = CGenericFactory::CreateTensorContainer(element);

	// Initialize the dimensions of C.
	C.resize(c, n2);
	
	// Compile-time tensor-product implementation.
	tensor->ResidualVolume(c, A1.data(), A2.data(), A3.data(), C.data());
	NTest_LA::MatrixErrorNormLinf(C, Bt, tol);
}


