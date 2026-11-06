#include "element/CStandardElement.hpp"
#include "common/AS3Constants.hpp"
#include "data/CMatrixAS3.hpp"
#include "quadrature/NLegendreGaussLobatto.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>



//-----------------------------------------------------------------------------------
// CStandardElement member functions.
//-----------------------------------------------------------------------------------


CStandardElement::CStandardElement
(
 ETypeNodalDistribution typeSolutionDistribution,
 std::size_t            nPolySolution,
 std::size_t            nBonusQuadratureFactor,
 ETypeQuadrature        typeQuadrature,
 ETypeDGMethod          typeDGMethod
)
	:
		mTypeDOFsSol(typeSolutionDistribution), mNPolySol(nPolySolution)
 /*
	* Constructor that initializes reference nodes, quadrature and basis operators.
	*/
{
	// Deduce the number of solution points in 1D and 2D.
	mNSol1D = mNPolySol + 1;
	mNSol2D = mNSol1D * mNSol1D;

	// Determine the required quadrature exactness for the configured DG method.
	std::size_t nPolyIntegrand = 0;
	switch( typeDGMethod )
	{
		case(ETypeDGMethod::DG_SEM):
		{
			if( nBonusQuadratureFactor   != 0 ||
			    typeSolutionDistribution != ETypeNodalDistribution::LGL ||
			    typeQuadrature           != ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO )
			{
				ERROR("DG-SEM requires zero bonus quadrature, LGL solution nodes and LGL quadrature.");
			}

			// Request LGL exactness through degree 2*p - 1 to obtain p + 1 points
			// collocated with the solution nodes. Degree-2*p basis products are underintegrated.
			nPolyIntegrand = 2 * nPolySolution - 1;
			break;
		}

		case(ETypeDGMethod::NODAL_DG):
		{
			// Products of two degree-p reference basis functions have degree 2*p.
			// The bonus adds nb*p as over-integration. Geometry and nonlinear fluxes
			// can require higher quadrature exactness.
			nPolyIntegrand = (2 + nBonusQuadratureFactor) * nPolySolution;
			break;
		}

		default: ERROR("Unknown DG method.");
	}


	// Initialize the integration rule and extract its reference nodes.
	mIntegrationRuleContainer = std::make_unique<CIntegrationRule>(typeQuadrature, nPolyIntegrand);
	const auto &rInt1D = mIntegrationRuleContainer->GetrInt1D();

	// Compute the solution DOFs in 1D and in parametric space.
	ComputeLocationDOFs1D(mTypeDOFsSol, mNSol1D, mRSol1D);

	// Collocation requires matching nodal distributions as well as matching point counts.
	mIsCollocated = mTypeDOFsSol == mIntegrationRuleContainer->GetTypeNodalDistribution() &&
	                mNSol1D == mIntegrationRuleContainer->GetnInt1D();



	// Initialize the Vandermonde matrices in 1D. Note, the bases are the orthonormal Legendre polynomials.
	NPolynomialUtility::OrthonormalLegendreBasis1D( mNSol1D, mRSol1D, mVandermondeSol1D ); // At solution    points.
	NPolynomialUtility::OrthonormalLegendreBasis1D( mNSol1D,  rInt1D, mVandermondeInt1D ); // At integration points.

	// Initialize the derivative of the Vandermonde matrices in 1D.
	NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D( mNSol1D, mRSol1D, mDerVandermondeSol1D ); // Solution    points.
	NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D( mNSol1D,  rInt1D, mDerVandermondeInt1D ); // Integration points.

	// Initialize the Lagrange polynomials in 1D.
	NPolynomialUtility::LagrangeBasis1D( mVandermondeInt1D, mVandermondeSol1D, mLagrangeInt1D );

	// Initialize the derivative of the Lagrange polynomials in 1D.
	NPolynomialUtility::DerivativeLagrangeBasis1D( mDerVandermondeSol1D, mVandermondeSol1D, mDerLagrangeSol1D );
	NPolynomialUtility::DerivativeLagrangeBasis1D( mDerVandermondeInt1D, mVandermondeSol1D, mDerLagrangeInt1D );

	// Initialize the derivative of the Lagrange polynomials at the solution min and max faces in 1D.  
	NPolynomialUtility::DerivativeLagrangeBasisFace1D( mDerLagrangeSol1D, mDerLagrangeMinFace1D, mDerLagrangeMaxFace1D );

	// Compute the transpose of the Lagrange polynomial at integration points in 1D.
	NLinearAlgebra::TransposeAS3Matrix(mLagrangeInt1D, mLagrangeInt1DTrans);
	// Compute the transpose of the derivative Lagrange polynomial at integration points in 1D.
	NLinearAlgebra::TransposeAS3Matrix(mDerLagrangeInt1D, mDerLagrangeInt1DTrans);
	// Compute the transpose of the derivative Lagrange polynomial at solution points in 1D.
	NLinearAlgebra::TransposeAS3Matrix(mDerLagrangeSol1D, mDerLagrangeSol1DTrans);

}

//-----------------------------------------------------------------------------------

CStandardElement::~CStandardElement
(
 void
)
 /*
	* Destructor for CStandardElement.
	*/
{

}

//-----------------------------------------------------------------------------------

void CStandardElement::ComputeLocationDOFs1D
(
 ETypeNodalDistribution typeDOFs,
 std::size_t            nDOFs1D,
 CMatrixAS3<as3double>  &rDOFs1D
)
 /*
	* Function that computes the requested reference nodes on [-1,+1].
	*/
{
	if( nDOFs1D == 0 ) ERROR("At least one solution point is required.");

	// Allocate memory for the nodes.
	rDOFs1D.resize(nDOFs1D);

	// Determine what type of distribution the DOFs obey.
	switch( typeDOFs )
	{
		// Equidistant distribution.
		case(ETypeNodalDistribution::EQD):
		{
			// A single equidistant point is the midpoint; avoid division by zero.
			if( nDOFs1D == 1 ) { rDOFs1D[0] = 0.0; break; }

			// Step size.
			const as3double dh = 2.0/( (as3double) nDOFs1D - 1.0 );

			// Starting position.
			const as3double r0 = -1.0;

			for(std::size_t i=0; i<nDOFs1D; i++)
				rDOFs1D[i] = (as3double) i*dh + r0;

			// Hard-enforce last node.
			rDOFs1D[nDOFs1D-1] = 1.0;

			break;
		}

		// Legendre-Gauss-Lobatto distribution.
		case(ETypeNodalDistribution::LGL):
		{
			// Preserve the midpoint convention for a constant solution.
			if( nDOFs1D == 1 ) { rDOFs1D[0] = 0.0; break; }

			if( NLegendreGaussLobatto::IsAvailable(nDOFs1D) )
			{
				const auto points = NLegendreGaussLobatto::GetPoints(nDOFs1D);
				for(std::size_t i=0; i<nDOFs1D; i++) rDOFs1D[i] = points[i];
			}
			else
			{
				// Request exactly nDOFs1D LGL points. The integration rule owns
				// numerical generation and its warning for untabulated counts.
				const CIntegrationRule rule(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 2 * (nDOFs1D - 2));
				const auto &points = rule.GetrInt1D();
				for(std::size_t i=0; i<nDOFs1D; i++) rDOFs1D[i] = points[i];
			}
			break;
		}

		// Gauss-Legendre solution nodes require face interpolation kernels.
		case(ETypeNodalDistribution::GL):
		{
			ERROR("Gauss-Legendre solution nodes are not supported by the current face kernels.");
			break;
		}

		default: ERROR("Unknown type of DOF distribution!");
	}

	// Make sure the DOFs are correct.
	as3double tmp = 0.0;
	for(std::size_t i=0; i<nDOFs1D; i++) tmp += rDOFs1D[i];
	const as3double roundoff  = static_cast<as3double>(nDOFs1D)*C_EPS;
	const as3double tolerance = std::max(C_ABS_TOL, roundoff);
	if( !(std::abs(tmp) <= tolerance) ) ERROR("Precision errors in nodal point computation.");
}


