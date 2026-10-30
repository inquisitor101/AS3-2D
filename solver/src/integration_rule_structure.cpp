#include "integration_rule_structure.hpp"

#include "error_structure.hpp"
#include "quadrature_structure.hpp"
#include "quadrature/gauss_legendre.hpp"
#include "quadrature/legendre_gauss_lobatto.hpp"

#include <string>


//-----------------------------------------------------------------------------------
// CIntegrationRule member functions.
//-----------------------------------------------------------------------------------


CIntegrationRule::CIntegrationRule
(
 ETypeQuadrature typeQuadrature,
 std::size_t     nPolyIntegrand
)
	:
		mTypeQuadrature(typeQuadrature)
 /*
	* Constructor that initializes one-dimensional quadrature data.
	*/
{
	switch( mTypeQuadrature )
	{
		case(ETypeQuadrature::GAUSS_LEGENDRE):
		{
			InitializeGaussLegendre(nPolyIntegrand);
			break;
		}

		case(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO):
		{
			InitializeLegendreGaussLobatto(nPolyIntegrand);
			break;
		}

		default: ERROR("Unknown quadrature rule.");
	}

  // Extract the number of integration points in 1D.
  const std::size_t nInt1D = GetnInt1D();

	// Ensure the generated data have the requested dimensions.
  if( nInt1D == 0 || nInt1D != mRInt1D.size() || nInt1D != mWInt1D.size() )
  {
		ERROR("Incorrect quadrature data dimensions.");
  }


  // Deduce the number of points in 2D.
  const std::size_t nInt2D = nInt1D * nInt1D;

	// Allocate the memory for the integration weights in 2D.
	mWInt2D.resize(nInt2D);

	// Form the integration weights in 2D, which are tensor-products of 1D vectors.
	std::size_t s = 0;
	for(std::size_t j=0; j<nInt1D; j++)
  {
		for(std::size_t i=0; i<nInt1D; i++)
    {
			mWInt2D[s++] = mWInt1D[i] * mWInt1D[j];
    }
  }
}

//-----------------------------------------------------------------------------------

std::size_t CIntegrationRule::GetnPolyExact1D
(
 void
) const
 /*
	* Function that returns the polynomial degree of exactness.
	*/
{
	const std::size_t nInt1D = GetnInt1D();

	switch( mTypeQuadrature )
	{
		case(ETypeQuadrature::GAUSS_LEGENDRE):         return 2*nInt1D - 1;
		case(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO): return 2*nInt1D - 3;

		default: ERROR("Unknown quadrature rule.");
	}

	// To avoid a compiler warning.
	return 0;
}

//-----------------------------------------------------------------------------------

std::size_t CIntegrationRule::ComputePointCount1D
(
 ETypeQuadrature typeQuadrature,
 std::size_t     nPolyIntegrand
)
 /*
	* Function that selects the point count without constructing a quadrature rule.
	*/
{
	switch( typeQuadrature )
	{
		case(ETypeQuadrature::GAUSS_LEGENDRE):         return nPolyIntegrand/2 + 1;
		case(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO): return nPolyIntegrand/2 + 2;

		default: ERROR("Unknown quadrature rule.");
	}
}

//-----------------------------------------------------------------------------------

void CIntegrationRule::InitializeGaussLegendre
(
 std::size_t nPolyIntegrand
)
 /*
	* Function that initializes Gauss-Legendre points and weights
	* on [-1,+1], which excludes both endpoints.
	*/
{
	// Select the point count required for the requested polynomial degree.
	const std::size_t nInt1D = ComputePointCount1D(ETypeQuadrature::GAUSS_LEGENDRE, nPolyIntegrand);

  // Consistency check.
  if( nInt1D < 1 )
  {
    ERROR("Gauss-Legendre quadrature requires at least one point.");
  }

	// Allocate memory for the integration points and weights in 1D.
	mRInt1D.resize(nInt1D);
	mWInt1D.resize(nInt1D);

	// Use the precomputed reference rule when this point count is tabulated.
	if( NGaussLegendre::IsAvailable(nInt1D) )
	{
		const auto points  = NGaussLegendre::GetPoints(nInt1D);
		const auto weights = NGaussLegendre::GetWeights(nInt1D);
		for(std::size_t i=0; i<nInt1D; i++)
		{
			mRInt1D[i] = points[i];
			mWInt1D[i] = weights[i];
		}
		return;
	}

  // Issue a warning, in case the requested rules are not precompiled.
	WARNING("Gauss-Legendre quadrature with " + std::to_string(nInt1D) +
	        " points is not tabulated. Generating with CGaussJacobiQuadrature; "
	        "this point count is not guaranteed to be supported.");

	// Allocate memory for temporary vectors, since they are double-precision only in 
	// the quadrature structure. NOTE, CGaussJacobiQuadrature must be refactored.
	as3vector1d<double> rtmp(nInt1D);
	as3vector1d<double> wtmp(nInt1D);

	// Compute the quadrature points and weights in 1D.
	CGaussJacobiQuadrature::GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rtmp, wtmp);

	// Copy the temporary data to their correct vector, in the proper precision.
	for(std::size_t i=0; i<nInt1D; i++)
	{
		mRInt1D[i] = static_cast<as3double>( rtmp[i] );
		mWInt1D[i] = static_cast<as3double>( wtmp[i] );
	}
}

//-----------------------------------------------------------------------------------

void CIntegrationRule::InitializeLegendreGaussLobatto
(
 std::size_t nPolyIntegrand
)
 /*
	* Function that initializes Legendre-Gauss-Lobatto points and weights
	* on [-1,+1], including both endpoints.
	*/
{
	// Select the point count required for the requested polynomial degree.
	const std::size_t nInt1D = ComputePointCount1D(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, nPolyIntegrand);

  // Consistency check.
	if( nInt1D < 2 )
  {
		ERROR("Legendre-Gauss-Lobatto quadrature requires at least two points.");
  }

	// Allocate the complete quadrature data.
	mRInt1D.resize(nInt1D);
	mWInt1D.resize(nInt1D);

	// Use the precomputed reference rule when this point count is tabulated.
	if( NLegendreGaussLobatto::IsAvailable(nInt1D) )
	{
		const auto points  = NLegendreGaussLobatto::GetPoints(nInt1D);
		const auto weights = NLegendreGaussLobatto::GetWeights(nInt1D);
		for(std::size_t i=0; i<nInt1D; i++)
		{
			mRInt1D[i] = points[i];
			mWInt1D[i] = weights[i];
		}
		return;
	}

  // Issue a warning, in case the requested rules are not precompiled.
	WARNING("Legendre-Gauss-Lobatto quadrature with " + std::to_string(nInt1D) +
	        " points is not tabulated. Generating with CGaussJacobiQuadrature; "
	        "this point count is not guaranteed to be supported.");

	// Compute the endpoint weight in double precision.
	const double n    = static_cast<double>(nInt1D);
	const double wEnd = 2.0/(n*(n-1.0));

	// Initialize the two endpoints.
	mRInt1D[0]        = -C_ONE;
	mRInt1D[nInt1D-1] =  C_ONE;

	mWInt1D[0]        = static_cast<as3double>(wEnd);
	mWInt1D[nInt1D-1] = static_cast<as3double>(wEnd);

	// The two-point rule contains only the endpoints.
	if( nInt1D == 2 ) return;

	// Number of interior quadrature points.
	const std::size_t nInterior = nInt1D - 2;

	// The existing quadrature generator operates in double precision.
	as3vector1d<double> rtmp(nInterior);
	as3vector1d<double> wtmp(nInterior);

	// Generate interior LGL points using Gauss-Jacobi quadrature
	// with weight function (1-r)*(1+r) = 1-r*r.
  // Note, CGaussJacobiQuadrature must be refactored.
	CGaussJacobiQuadrature::GetQuadraturePoints(1.0, 1.0, -1.0, 1.0, rtmp, wtmp);

	// Convert the weighted Jacobi rule to the unweighted LGL rule.
	for(std::size_t i=0; i<nInterior; i++)
	{
		const double r = rtmp[i];
		const double w = wtmp[i]/(1.0-r*r);

		mRInt1D[i+1] = static_cast<as3double>(r);
		mWInt1D[i+1] = static_cast<as3double>(w);
	}
}
