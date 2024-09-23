#include "polynomial_structure.hpp"


//-----------------------------------------------------------------------------------
// NPolynomialUtility namespace functions.
//-----------------------------------------------------------------------------------


unsigned short NPolynomialUtility::IntegrationRule1D
(
 unsigned short npoly
)
 /*
	* Function that estimates the number of integration points in 1D.
	*/
{
	// The (over-)integration rule specified.
	const unsigned short nd = 3; 	
	
	// Cast the value of the number of integration points, then return it.
	return static_cast<unsigned short>(nd*npoly/2 + 1);
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::LagrangeBasis1D
(
 CMatrixAS3<as3double> &V1,
 CMatrixAS3<as3double> &V2,
 CMatrixAS3<as3double> &L1D
)
 /*
	* Function that computes the the Lagrange basis functions in 1D via 2 Vandermonde 
	* matrices: L = V1*inv(V2).
	*/
{
	// Determine the number of rows and columns in the Lagrange polynomial.
	const size_t nRow = V1.row();
	const size_t nCol = V1.col();

	// Consistency check for matrix dimensions. Note, V2 must be a square matrix.
	if( (nCol != V2.row()) || (nCol != V2.col()) ) ERROR("Matrix dimensions are incorrect.");

	// Temporary data that stores the inverse of matrix V2.
	CMatrixAS3<as3double> V2inv = V2;

	// Compute the inverse of the Vandermonde matrix V2.
	NLinearAlgebra::InverseMatrix( V2inv ); 

	// Compute the Lagrange basis function.
	NLinearAlgebra::MatrixMatrixMult( V1, V2inv, L1D ); 

	// Ensure the correct dimensions of the resulting matrix.
	if( (L1D.row() != nRow) || (L1D.col() != nCol) ) ERROR("Result of Lagrange basis is incorrect.");

	// Check that the Lagrange basis functions are correct, by ensuring the sum of each row is unity.
	as3double one = 1.0;	
	for(size_t i=0; i<nRow; i++)
	{
		as3double tmp = 0.0;
		for(size_t j=0; j<nCol; j++)  tmp += L1D(i,j);
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
	* Function that computes the Lagrange differentiation matrix in 1D via 2 Vandermonde 
	* matrices: dL = dV1*inv(V2).
	*/
{
	// Determine the number of rows and columns in the derivative of the Lagrange polynomial.
	const size_t nRow = dV1.row();
	const size_t nCol = dV1.col();

	// Consistency check for matrix dimensions. Note, V2 must be a square matrix.
	if( (nCol != V2.row()) || (nCol != V2.col()) ) ERROR("Matrix dimensions are incorrect.");

	// Temporary data that stores the inverse of matrix V2.
	CMatrixAS3<as3double> V2inv = V2;

	// Compute the inverse of the Vandermonde matrix V2.
	NLinearAlgebra::InverseMatrix( V2inv ); 

	// Compute the derivative of the Lagrange basis function.
	NLinearAlgebra::MatrixMatrixMult( dV1, V2inv, dL1D ); 

	// Ensure the correct dimensions of the resulting matrix.
	if( (dL1D.row() != nRow) || (dL1D.col() != nCol) ) ERROR("Result of differentiation matrix is incorrect.");

	// Check that the differentiation matrix is correct, by ensuring the sum of each row is zero.
	for(size_t i=0; i<nRow; i++)
	{
		as3double tmp = 0.0;
		for(size_t j=0; j<nCol; j++) tmp += dL1D(i,j);
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
	* Function that computes the derivative of the Lagrange basis function 
	* at the min (r = -1) and max (r = +1) faces in 1D.
	*/
{
	// Ensure the input matrix is indeed square.
	if( dLSol1D.row() != dLSol1D.col() ) ERROR("Input differentiation matrix must be square.");

	// Explicitly deduce the dimension of the solution points.
	const size_t n = dLSol1D.row();

	// Allocate memory for both faces.
	dLmin.resize( n );
	dLmax.resize( n );

	// Compute the derivative coefficients at the faces.
	for(size_t i=0; i<n; i++)
	{
		dLmin[i] = dLSol1D(0,  i);
		dLmax[i] = dLSol1D(n-1,i); 
	}

	// Check that the differentiation coefficients are correct, by ensuring the sum of each is zero.
	as3double min = 0.0, max = 0.0;
	for(size_t i=0; i<n; i++)
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
 const unsigned int           nCol,
 const CMatrixAS3<as3double> &rDOFs1D,
 CMatrixAS3<as3double>       &rBasis1D
)
 /*
	* Function that computes the Legendre basis functions in 1D.
	*/
{
	// Determine the number of rows.
	const size_t nRow = rDOFs1D.size();

	// Allocate memory for the output basis.
	rBasis1D.resize(nRow, nCol);

	for(unsigned int iRow=0; iRow<nRow; iRow++)
		for(unsigned int iCol=0; iCol<nCol; iCol++)
			rBasis1D(iRow,iCol) = 
				NormalizedLegendrePolynomial(static_cast<int>(iCol), rDOFs1D[iRow]);
}

//-----------------------------------------------------------------------------------

void NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D
(
 const unsigned int           nCol,
 const CMatrixAS3<as3double> &rDOFs1D,
 CMatrixAS3<as3double>       &rDerBasis1D
)
 /*
	* Function that computes the derivative of the Legendre basis functions in 1D.
	*/
{
	// Determine the number of rows.
	const size_t nRow = rDOFs1D.size();

	// Allocate memory for the output basis.
	rDerBasis1D.resize(nRow, nCol);

	for(unsigned int iRow=0; iRow<nRow; iRow++)
		for(unsigned int iCol=0; iCol<nCol; iCol++)
			rDerBasis1D(iRow,iCol) = 
				DerivativeNormalizedLegendrePolynomial(static_cast<int>(iCol), rDOFs1D[iRow]);
}

//-----------------------------------------------------------------------------------

as3double NPolynomialUtility::NormalizedLegendrePolynomial
(
 const int       n,
 const as3double r
)
 /*
	* Function that computes the coefficient of the normalized Legendre polynomial at a point.
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
 const int       n,
 const as3double r
)
 /*
	* Function that computes the coefficient of the derivative of the normalized 
	* Legendre polynomial at a point.
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
 const int       n, 
 const int       alpha,
 const int       beta,
 const as3double r
)
 /*
	* Function that computes the value of the derivative of the normalized Jacobi 
	* polynomial with coefficients alpha and beta of order n for the given value of r.
	*
	* NOTE, this function is taken from Edwin van der Weide's implementatoin in VCP3D.
	*/
{
  // Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
  // zero, because the polynomial itself is constant.
  as3double grad;
  if(n == 0) grad = 0.0;
  else
  {
    const as3double tmp = n*(n+alpha+beta+1.0);
    grad = std::sqrt(tmp)*NormalizedJacobiPolynomial(n-1, alpha+1, beta+1, r);
  }

  // Return the gradient.
  return grad;
}

//-----------------------------------------------------------------------------------

as3double NPolynomialUtility::NormalizedJacobiPolynomial
(
 const int       n, 
 const int       alpha,
 const int       beta,
 const as3double r
)
 /*
	* Function that computes the value of the normalized Jacobi polynomial with 
	* coefficients alpha and beta of order n for the given value of r.
	*
	* NOTE, this function is taken from Edwin van der Weide's implementatoin in VCP3D.
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
    for(int i=2; i<=n; ++i)
    {
      // Compute the coefficients a for i and i-1 and the coefficient bi.
      int j = i-1;
      as3double tmp  = 2*j + apb;
      as3double aim1 = 2.0*std::sqrt(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-1.0)*(tmp+1.0)))/tmp;

      as3double bi = b2ma2/(tmp*(tmp+2.0));

      tmp = 2*i + apb;
      as3double ai = 2.0*std::sqrt(i*(i+apb)*(i+alpha)*(i+beta)/((tmp-1.0)*(tmp+1.0)))/tmp;

      // Compute the new value of Pn and make sure to store Pnm1 correctly.
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((r-bi)*Pn - aim1*tmp)/ai;
    }
  }

  // Return Pn.
  return Pn;
}


