#include "polynomial_structure.hpp"

#include <cmath>
#include <cstddef>


//-----------------------------------------------------------------------------------
// NPolynomialUtility namespace functions.
//-----------------------------------------------------------------------------------


std::size_t NPolynomialUtility::IntegrationRule1D
(
 std::size_t npoly
)
 /*
	* Function that selects the one-dimensional over-integration point count.
	*/
{
	// The (over-)integration rule specified.
	const std::size_t nd = 3; 	
	
	// Return the number of integration points without narrowing it.
	return nd*npoly/2 + 1;
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::LagrangeBasis1D
(
 CMatrixAS3<as3double> &V1,
 CMatrixAS3<as3double> &V2,
 CMatrixAS3<as3double> &L1D
)
 /*
	* Function that forms the Lagrange interpolation matrix as V1*inv(V2).
	*/
{
	// Determine the number of rows and columns in the Lagrange polynomial.
	const std::size_t nRow = V1.row();
	const std::size_t nCol = V1.col();

	// Consistency check for matrix dimensions. Note, V2 must be a square matrix.
	if( (nCol != V2.row()) || (nCol != V2.col()) ) ERROR("Matrix dimensions are incorrect.");

	// Temporary data that stores the inverse of matrix V2.
	CMatrixAS3<as3double> V2inv = V2;

	// Compute the inverse of the Vandermonde matrix V2.
	NLinearAlgebra::InverseMatrix( V2inv ); 

	// Compute the Lagrange basis function.
	L1D = NLinearAlgebra::MatrixMatrixMult(V1, V2inv); 

	// Ensure the correct dimensions of the resulting matrix.
	if( (L1D.row() != nRow) || (L1D.col() != nCol) ) ERROR("Result of Lagrange basis is incorrect.");

	// Check that the Lagrange basis functions are correct, by ensuring the sum of each row is unity.
	as3double one = 1.0;	
	for(std::size_t i=0; i<nRow; i++)
	{
		as3double tmp = 0.0;
		for(std::size_t j=0; j<nCol; j++)  tmp += L1D(i,j);
		if( fabs(tmp-one) > 1.0e-10 ) ERROR("Results of the Lagrange interpolation function are incorrect.");
	}
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::DerivativeLagrangeBasis1D
(
 CMatrixAS3<as3double> &dV1,
 CMatrixAS3<as3double> &V2,
 CMatrixAS3<as3double> &dL1D
)
 /*
	* Function that forms the Lagrange differentiation matrix as dV1*inv(V2).
	*/
{
	// Determine the number of rows and columns in the derivative of the Lagrange polynomial.
	const std::size_t nRow = dV1.row();
	const std::size_t nCol = dV1.col();

	// Consistency check for matrix dimensions. Note, V2 must be a square matrix.
	if( (nCol != V2.row()) || (nCol != V2.col()) ) ERROR("Matrix dimensions are incorrect.");

	// Temporary data that stores the inverse of matrix V2.
	CMatrixAS3<as3double> V2inv = V2;

	// Compute the inverse of the Vandermonde matrix V2.
	NLinearAlgebra::InverseMatrix( V2inv ); 

	// Compute the derivative of the Lagrange basis function.
	dL1D = NLinearAlgebra::MatrixMatrixMult(dV1, V2inv); 

	// Ensure the correct dimensions of the resulting matrix.
	if( (dL1D.row() != nRow) || (dL1D.col() != nCol) ) ERROR("Result of differentiation matrix is incorrect.");

	// Check that the differentiation matrix is correct, by ensuring the sum of each row is zero.
	for(std::size_t i=0; i<nRow; i++)
	{
		as3double tmp = 0.0;
		for(std::size_t j=0; j<nCol; j++) tmp += dL1D(i,j);
		if( fabs(tmp) > 1.0e-10 )    ERROR("Results of the Lagrange differentiation matrix are incorrect.");
	}
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::DerivativeLagrangeBasisFace1D
(
 const CMatrixAS3<as3double> &dLSol1D,
 CMatrixAS3<as3double>       &dLmin,
 CMatrixAS3<as3double>       &dLmax
)
 /*
	* Function that extracts basis derivatives at the minimum and maximum reference endpoints.
	*/
{
	// Ensure the input matrix is indeed square.
	if( dLSol1D.row() != dLSol1D.col() ) ERROR("Input differentiation matrix must be square.");

	// Explicitly deduce the dimension of the solution points.
	const std::size_t n = dLSol1D.row();

	// Allocate memory for both faces.
	dLmin.resize( n );
	dLmax.resize( n );

	// Compute the derivative coefficients at the faces.
	for(std::size_t i=0; i<n; i++)
	{
		dLmin[i] = dLSol1D(0,  i);
		dLmax[i] = dLSol1D(n-1,i); 
	}

	// Check that the differentiation coefficients are correct, by ensuring the sum of each is zero.
	as3double min = 0.0, max = 0.0;
	for(std::size_t i=0; i<n; i++)
	{
		min += dLmin[i];
		max += dLmax[i];
	}

	// Relative error tolerance.
	const as3double tol = 1.0e-10;
	if( (fabs(min) > tol) || (fabs(max) > tol) ) 
		ERROR("Results of the Lagrange differentiation at face are incorrect.");
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::OrthonormalLegendreBasis1D
(
 const std::size_t           nCol,
 const CMatrixAS3<as3double> &rDOFs1D,
 CMatrixAS3<as3double>       &rBasis1D
)
 /*
	* Function that evaluates orthonormal Legendre basis functions at the supplied nodes.
	*/
{
	// Determine the number of rows.
	const std::size_t nRow = rDOFs1D.size();

	// Allocate memory for the output basis.
	rBasis1D.resize(nRow, nCol);

	for(std::size_t iRow=0; iRow<nRow; iRow++)
		for(std::size_t iCol=0; iCol<nCol; iCol++)
			rBasis1D(iRow,iCol) = 
				NormalizedLegendrePolynomial(iCol, rDOFs1D[iRow]);
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D
(
 const std::size_t           nCol,
 const CMatrixAS3<as3double> &rDOFs1D,
 CMatrixAS3<as3double>       &rDerBasis1D
)
 /*
	* Function that evaluates derivatives of the orthonormal Legendre basis at the supplied nodes.
	*/
{
	// Determine the number of rows.
	const std::size_t nRow = rDOFs1D.size();

	// Allocate memory for the output basis.
	rDerBasis1D.resize(nRow, nCol);

	for(std::size_t iRow=0; iRow<nRow; iRow++)
		for(std::size_t iCol=0; iCol<nCol; iCol++)
			rDerBasis1D(iRow,iCol) = 
				DerivativeNormalizedLegendrePolynomial(iCol, rDOFs1D[iRow]);
}

//-----------------------------------------------------------------------------------

as3double NPolynomialUtility::NormalizedLegendrePolynomial
(
 const std::size_t       n,
 const as3double r
)
 /*
	* Function that evaluates the normalized Legendre polynomial of the requested order.
	*/
{
	// Ensure the evaluation range is between: [-1,+1].
	const as3double rmin = (as3double) -1.0;
	const as3double rmax = (as3double)  1.0;
	if( (r < rmin) || (r > rmax) ) ERROR("Evaluation point outside interval: [-1,+1].");

	// Recall, the normalized Legendre polynomials correspond to the special case of the 
	// normalized Jacobi polynomials with: alpha = beta = 0.
	return NormalizedJacobiPolynomial(n, 0, 0, r); 	
}

//-----------------------------------------------------------------------------------

as3double NPolynomialUtility::DerivativeNormalizedLegendrePolynomial
(
 const std::size_t       n,
 const as3double r
)
 /*
	* Function that evaluates the derivative of a normalized Legendre polynomial.
	*/
{
	// Ensure the evaluation range is between: [-1,+1].
	const as3double rmin = (as3double) -1.0;
	const as3double rmax = (as3double)  1.0;
	if( (r < rmin) || (r > rmax) ) ERROR("Evaluation point outside interval: [-1,+1].");

	// Recall, the normalized Legendre polynomials correspond to the special case of the 
	// normalized Jacobi polynomials with: alpha = beta = 0.
	return DerivativeNormalizedJacobiPolynomial(n, 0, 0, r); 	
}

//-----------------------------------------------------------------------------------

as3double NPolynomialUtility::DerivativeNormalizedJacobiPolynomial
(
 const std::size_t       n, 
 const int       alpha,
 const int       beta,
 const as3double r
)
 /*
	* Function that evaluates the derivative of a normalized Jacobi polynomial.
	*
	* NOTE, adapted from Edwin van der Weide's implementation in VCP3D.
	*/
{
  // Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
  // zero, because the polynomial itself is constant.
  as3double grad;
  if(n == 0) grad = 0.0;
  else
  {
    const as3double degree = static_cast<as3double>(n);
    const as3double tmp = degree*(degree+alpha+beta+1.0);
    grad = std::sqrt(tmp)*NormalizedJacobiPolynomial(n-1, alpha+1, beta+1, r);
  }

  // Return the gradient.
  return grad;
}

//-----------------------------------------------------------------------------------

as3double NPolynomialUtility::NormalizedJacobiPolynomial
(
 const std::size_t       n, 
 const int       alpha,
 const int       beta,
 const as3double r
)
 /*
	* Function that evaluates a normalized Jacobi polynomial at the supplied point.
	*
	* NOTE, adapted from Edwin van der Weide's implementation in VCP3D.
	*/
{
	// Ensure the evaluation range is between: [-1,+1].
	const as3double rmin = (as3double) -1.0;
	const as3double rmax = (as3double)  1.0;
	if( (r < rmin) || (r > rmax) ) ERROR("Evaluation point outside interval: [-1,+1].");

  // Some abbreviations.
  const as3double ap1   = (as3double) (alpha + 1);
  const as3double bp1   = (as3double) (beta  + 1);
  const as3double apb   = (as3double) (alpha + beta);
  const as3double apbp1 = (as3double) (apb + 1);
  const as3double apbp2 = (as3double) (apb + 2);
  const as3double apbp3 = (as3double) (apb + 3);
  const as3double b2ma2 = (as3double) (beta*beta - alpha*alpha);

  // Determine the terms, which involves the gamma function. As the
  // arguments are integers, this term can be computed easily, because
  // Gamma(n+1) = n!.
  as3double Gamap1 = 1.0, Gambp1 = 1.0, Gamapbp2 = 1.0;
  for(int i=2; i<=alpha; ++i)          Gamap1   *= i;
  for(int i=2; i<=beta; ++i)           Gambp1   *= i;
  for(int i=2; i<=(alpha+beta+1); ++i) Gamapbp2 *= i;

  // Initialize the normalized polynomials.
  as3double Pnm1 = std::sqrt(pow(0.5,apbp1)*Gamapbp2/(Gamap1*Gambp1));
  as3double Pn   = 0.5*Pnm1*(apbp2*r + alpha - beta)*std::sqrt(apbp3/(ap1*bp1));

  // Take care of the special situation of n == 0.
  if(n == 0) Pn = Pnm1;
  else
  {
    // The value of the normalized Jacobi polynomial must be obtained
    // via recursion.
    for(std::size_t i=2; i<=n; ++i)
    {
      // Compute the coefficients a for i and i-1 and the coefficient bi.
      const as3double j = static_cast<as3double>(i-1);
      const as3double degree = static_cast<as3double>(i);
      as3double tmp  = 2*j + apb;
      as3double aim1 = 2.0*std::sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))/tmp;

      as3double bi = b2ma2/(tmp*(tmp+2.0));

      tmp = 2*degree + apb;
      as3double ai = 2.0*std::sqrt(degree*(degree+apb)*(degree+alpha)*(degree+beta)/((tmp-1.0)*(tmp+1.0)))/tmp;

      // Compute the new value of Pn and make sure to store Pnm1 correctly.
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((r-bi)*Pn - aim1*tmp)/ai;
    }
  }

  // Return Pn.
  return Pn;
}


