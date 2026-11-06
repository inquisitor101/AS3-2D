#include "solver/CEESinglezoneSolver.hpp"
#include "data/CPoolMatrixAS3.hpp"

#include "compute_kernels/TensorProducts.hpp"
#include "CGenericFactory.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "element/CPhysicalElement.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "element/CStandardElement.hpp"

#include <cstddef>
#include <iostream>


//-----------------------------------------------------------------------------------
// CEESinglezoneSolver member functions.
//-----------------------------------------------------------------------------------


CEESinglezoneSolver::CEESinglezoneSolver
(
 const CConfig             *config_container,
 const CSinglezoneGeometry *singlezone_geometry_container
)
	:
		ISinglezoneSolver(config_container, singlezone_geometry_container, CEESinglezoneSolver::mNVar)
 /*
	* Constructor that initializes the Euler solver for one zone.
	*/
{

}

//-----------------------------------------------------------------------------------

CEESinglezoneSolver::~CEESinglezoneSolver
(
 void
)
 /*
	* Destructor for CEESinglezoneSolver.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::InitPhysicalElements
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that constructs physical element states and metrics for the zone.
	*/
{	
	// Report output.
	std::cout << "    physical elements.... "; 

	// Get a reference to the current zone.
	auto* zone = multizone_geometry_container->GetSinglezoneGeometry(mZoneID);

	// Allocate memory for the physical elements in this solver.
	mPhysicalElementContainer.resize( zone->GetnElem() );
	
	// Instantiate physical element objects. The reason this needs a factory is because we might have 
	// a specialized solver (e.g. EE), but with a buffer zone, such as a PML. Hence, the elements
	// require additional variables. For now, avoid an interface class, as this object is the most 
	// performance-reliant data class in the code. Perhaps it is better to keep a single class, but 
	// instantiate it differently, depending on the type of solver and buffer layer specification.
	for(std::size_t i=0; i<mPhysicalElementContainer.size(); i++)
	{
		mPhysicalElementContainer[i] = CGenericFactory::CreatePhysicalElement(config_container,
				                                                                  mStandardElementContainer.get(),
																																					mTensorProductContainer.get(),
				                                                                  zone->GetElementGeometry(i),
																																					mZoneID, mNVar);
	}

	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::ComputeVolumeResidual
(
 CSinglezoneGeometry       *singlezone_geometry_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElem
)
 /*
	* Function that overwrites one element residual with its volume contribution.
	*/
{
	// Extract the number of integration points in 2D.
	const std::size_t nInt2D = mStandardElementContainer->GetnInt2D();
	// Extract the quadrature integration weights on the standard element.
	auto&  wInt2D = mStandardElementContainer->GetwInt2D();

	// Borrow memory for the solution and its gradient.
	CWorkMatrixAS3<as3double> var    = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
	CWorkMatrixAS3<as3double> dVarDx = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
	CWorkMatrixAS3<as3double> dVarDy = workarray.GetWorkMatrixAS3(mNVar, nInt2D);

	// Extract the relevant physical element.
	auto& physical_element = mPhysicalElementContainer[iElem];

	// Extract the metrics at the volume solution points.
	auto& jac = physical_element->mMetricInt2D;
	// Reference to the current element solution.
	auto& sol = physical_element->mSol2D;
	// Reference to the current element residual.
	auto& res = physical_element->mRes2D;

	// Compute the solution and its (parametric) gradient at the volume integration points.
	mTensorProductContainer->CompileTimeVolume(sol.data(),
			                                       var.data(), dVarDx.data(), dVarDy.data());

	// Convert the gradient from parametric to Cartesian coordinates.
	physical_element->ConvertGradParamToCartVolInt(dVarDx, dVarDy);


	// Loop over the integration points and compute the weighted terms.
#pragma omp simd
  for(std::size_t l=0; l<nInt2D; l++)
	{
		// Compute the weighting factor: Jacobian multiplied by the integration weight.
		// Note, the residual is defined on the RHS: 
		//  dUDt = Res(vol) - Res(surf), thus the vol terms are +ve.
		const as3double weight = jac(0,l)*wInt2D[l];

		// Assemble the relevant (weighted) metrics in the x- and y-derivatives.
		const as3double wdrdx = weight*jac(1,l);
		const as3double wdrdy = weight*jac(2,l);
		const as3double wdsdx = weight*jac(3,l);
		const as3double wdsdy = weight*jac(4,l);

	  // Compute the primitive variables.
	  const as3double rho   = var(0,l);
	  const as3double ovrho = C_ONE/rho;
	  const as3double u     = ovrho*var(1,l);
	  const as3double v     = ovrho*var(2,l);
	  const as3double p     = C_GM1*( var(3,l) - C_HALF*(u*var(1,l) + v*var(2,l)) );

	  // Compute the inviscid flux in the x-direction. 
	  const as3double fx0 =     var(1,l);
	  const as3double fx1 = u*  var(1,l) + p;
	  const as3double fx2 = v*  var(1,l);
	  const as3double fx3 = u*( var(3,l) + p );

	  // Compute the inviscid flux in the y-direction. 
	  const as3double fy0 =     var(2,l);
	  const as3double fy1 = u*  var(2,l);
	  const as3double fy2 = v*  var(2,l) + p;
	  const as3double fy3 = v*( var(3,l) + p );

		// Compute the residual flux terms in the r-parametric direction.
		dVarDx(0,l) = fx0*wdrdx + fy0*wdrdy;
		dVarDx(1,l) = fx1*wdrdx + fy1*wdrdy;
		dVarDx(2,l) = fx2*wdrdx + fy2*wdrdy;
		dVarDx(3,l) = fx3*wdrdx + fy3*wdrdy;

		// Compute the residual flux terms in the s-parametric direction.
		dVarDy(0,l) = fx0*wdsdx + fy0*wdsdy;
		dVarDy(1,l) = fx1*wdsdx + fy1*wdsdy;
		dVarDy(2,l) = fx2*wdsdx + fy2*wdsdy;
		dVarDy(3,l) = fx3*wdsdx + fy3*wdsdy;
	}

	// Scatter the volume residual terms back to the actual residual. Note, this resets the residual! 
	mTensorProductContainer->CompileTimeResidualVolume(nullptr, 
                                                     dVarDx.data(), 
                                                     dVarDy.data(), 
                                                     res.data());
}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::ComputeSurfaceResidualIDir
(
 const CSinglezoneGeometry *singlezone_geometry_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElemL
)
 /*
	* Function that computes one same-zone I-face flux and fills both local face buffers.
	*/
{
  // Deduce the right element's index.
  const std::size_t iElemR = iElemL + 1;

	// Extract the number of integration points in 1D.
	std::size_t nInt1D = mStandardElementContainer->GetnInt1D();
	// Extract the quadrature integration weights in 1D on the standard element.
	auto&  wInt1D = mStandardElementContainer->GetwInt1D();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varL = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> varR = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, nInt1D);

	// Get a pointer to the respective left and right element, w.r.t. this face.
	auto& elemL = mPhysicalElementContainer[iElemL];
	auto& elemR = mPhysicalElementContainer[iElemR];

	// Extract the metrics at the integration points on the imax face of the left element.
	auto& metL = elemL->mMetricIntIMax1D;
	// Reference to the left element solution.
	auto& solL = elemL->mSol2D;
	// Reference to the right element residual (temporary).
	auto& resL = elemL->GetSurfaceResidual(EFaceLocation::IMAX);

	// Reset the values in the resL, before proceeding.
	for(std::size_t l=0; l<resL.size(); l++) resL[l] = C_ZERO;

	// Reference to the right element solution.
	auto& solR = elemR->mSol2D;
	// Reference to the right element residual.
	auto& resR = elemR->GetSurfaceResidual(EFaceLocation::IMIN);
  resR.reset();

	// Compute the solution on the integration nodes of the left element.
	mTensorProductContainer->CompileTimeSurfaceIMAX(solL.data(),
                                                  varL.data(), nullptr, nullptr); 

	// Compute the solution on the integration nodes of the right element.
	mTensorProductContainer->CompileTimeSurfaceIMIN(solR.data(),
		                                              varR.data(), nullptr, nullptr); 

	// Compute the flux state, weighted by the integration nodes and metrics. 
	// Notice, this is based on the left state, which is the outward-pointing
	// normal vector.
	mRiemannSolverContainer->ComputeFlux(wInt1D, metL, varL, varR, flux);

	// Compute the residual on the left  element, which is on the IMAX boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceIMAX(flux.data(), 
                                                          nullptr, nullptr, resL.data());

	// For local conservation, negate the flux, since it leaves the left element to enter the right element.
	for(std::size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

	// Compute the residual on the right element, which is on the IMIN boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceIMIN(flux.data(), 
                                                          nullptr, nullptr, resR.data());
}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::ComputeSurfaceResidualJDir
(
 const CSinglezoneGeometry *singlezone_geometry_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElemB
)
 /*
	* Function that computes one same-zone J-face flux and fills both local face buffers.
	*/
{
  // Extract the number of elements in the ith-direction in this zone.
	const std::size_t niElem = singlezone_geometry_container->GetniElem();	
  // Deduce the top element's index.
  const std::size_t iElemT = iElemB + niElem;

	// Extract the number of integration points in 1D.
	std::size_t nInt1D = mStandardElementContainer->GetnInt1D();
	// Extract the quadrature integration weights in 1D on the standard element.
	auto&  wInt1D = mStandardElementContainer->GetwInt1D();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varB = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> varT = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, nInt1D);

	// Get a pointer to the respective left and right element, w.r.t. this face.
	auto& elemB = mPhysicalElementContainer[iElemB];
	auto& elemT = mPhysicalElementContainer[iElemT];

	// Extract the metrics at the integration points on the imax face of the left element.
	auto& metB = elemB->mMetricIntJMax1D;
	// Reference to the left element solution.
	auto& solB = elemB->mSol2D;
	// Reference to the right element residual (temporary).
	auto& resB = elemB->GetSurfaceResidual(EFaceLocation::JMAX);

	// Reset the values in the resL, before proceeding.
	for(std::size_t l=0; l<resB.size(); l++) resB[l] = C_ZERO;

	// Reference to the right element solution.
	auto& solT = elemT->mSol2D;
	// Reference to the right element residual.
	auto& resT = elemT->GetSurfaceResidual(EFaceLocation::JMIN);
  resT.reset();

	// Compute the solution on the integration nodes of the left element.
	mTensorProductContainer->CompileTimeSurfaceJMAX(solB.data(),
		                                              varB.data(), nullptr, nullptr); 

	// Compute the solution on the integration nodes of the right element.
	mTensorProductContainer->CompileTimeSurfaceJMIN(solT.data(),
		                                              varT.data(), nullptr, nullptr); 

	// Compute the flux state, weighted by the integration nodes and metrics. 
	// Notice, this is based on the left state, which is the outward-pointing
	// normal vector.
	mRiemannSolverContainer->ComputeFlux(wInt1D, metB, varB, varT, flux);

	// Compute the residual on the left  element, which is on the IMAX boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceJMAX(flux.data(), 
                                                          nullptr, nullptr, resB.data());

	// For local conservation, negate the flux, since it leaves the left element to enter the right element.
	for(std::size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

	// Compute the residual on the right element, which is on the IMIN boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceJMIN(flux.data(), 
                                                          nullptr, nullptr, resT.data());
}
