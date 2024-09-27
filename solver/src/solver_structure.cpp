#include "solver_structure.hpp"


//-----------------------------------------------------------------------------------
// ISolver member functions.
//-----------------------------------------------------------------------------------


ISolver::ISolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
	:
		mZoneID(iZone)
 /*
	* Constructor for the interface solver class.
	*/
{
	// Instantiate a standard element object.
	mStandardElementContainer = CGenericFactory::CreateStandardElement(config_container, iZone);

	// Instantiate a tensor-product object.
	mTensorProductContainer = CGenericFactory::CreateTensorContainer(mStandardElementContainer.get());

	// Instantiate a Riemann solver object.
	mRiemannSolverContainer = CGenericFactory::CreateRiemannSolverContainer(config_container,
			                                                                    config_container->GetTypeRiemannSolver(iZone));

	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(iZone);

	// Ensure the zone matches this one.
	if( iZone != zone->GetZoneID() ) ERROR("Geometry zone does not match Solver zone.");
}

//-----------------------------------------------------------------------------------

ISolver::~ISolver
(
 void
)
 /*
	* Destructor, which cleans up after the interface solver class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CEESolver member functions.
//-----------------------------------------------------------------------------------


CEESolver::CEESolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
	:
		ISolver(config_container, geometry_container, iZone)
 /*
	* Constructor for the (non-linear) Euler equations class.
	*/
{

}

//-----------------------------------------------------------------------------------

CEESolver::~CEESolver
(
 void
)
 /*
	* Destructor, which cleans up after the Euler equations class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEESolver::InitPhysicalElements
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the physical elements. 
	*/
{	
	// Report output.
	std::cout << "    physical elements.... "; 

	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(mZoneID);

	// Allocate memory for the physical elements in this solver.
	mPhysicalElementContainer.resize( zone->GetnElem() );
	
	// Instantiate physical element objects. The reason this needs a factory is because we might have 
	// a specialized solver (e.g. EE), but with a buffer zone, such as a PML. Hence, the elements
	// require additional variables. For now, avoid an interface class, as this object is the most 
	// performance-reliant data class in the code. Perhaps it is better to keep a single class, but 
	// instantiate it differently, depending on the type of solver and buffer layer specification.
	for(size_t i=0; i<mPhysicalElementContainer.size(); i++)
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

void CEESolver::InitBoundaryConditions
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the boundary conditions. 
	*/
{
	// Report output.
	std::cout << "    boundary conditions.. ";

	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(mZoneID);

	// Index counters, used for determining the number of boundaries needed.
	size_t nimin = 0; size_t njmin = 0;
	size_t nimax = 0; size_t njmax = 0;

	// Initialize the boundary conditions, per each element face.
	
	// First, loop over the markers in this zone to determine how much memory we need
	// per boundary face container. Otherwise, we end up with capacity > size, which 
	// is not exactly efficient. 
	for( auto& marker: zone->GetMarker() )
	{
		// In case this marker is an interface, skip it.
		if( marker->GetTypeBC() == ETypeBC::INTERFACE ) continue;

		// Loop over each element on this marker and instantiate its boundary.
		for( auto& [index, face]: marker->GetElementFaces() )
		{
			// Make sure that the marker data is written in the correct order.
			static_assert( std::is_same_v<const unsigned int, decltype(index)> );
			static_assert( std::is_same_v<const EFaceElement, decltype(face )> );

			// Based on the face type, initialzie the boundary surface.
			switch( face )
			{

				// Check the IMIN boundary.
				case(EFaceElement::IMIN): 
				{
					// If the number of faces exceeds the expected, issue an error.
					if( nimin++ > zone->GetnyElem() ) ERROR("Incorrect number of IMIN faces."); break;
				}

				// Check the IMAX boundary.
				case(EFaceElement::IMAX): 
				{
					// If the number of faces exceeds the expected, issue an error.
					if( nimax++ > zone->GetnyElem() ) ERROR("Incorrect number of IMAX faces."); break;
				}

				// Check the JMIN boundary.
				case(EFaceElement::JMIN): 
				{
					// If the number of faces exceeds the expected, issue an error.
					if( njmin++ > zone->GetnxElem() ) ERROR("Incorrect number of JMIN faces."); break;
				}

				// Check the JMAX boundary.
				case(EFaceElement::JMAX): 
				{
					// If the number of faces exceeds the expected, issue an error.
					if( njmax++ > zone->GetnxElem() ) ERROR("Incorrect number of JMAX faces."); break;
				}

				// Something went wrong, issue an error.
				default: ERROR("Unknown face detected.");
			}
		}
	}

	// Reserve memory for the boundary faces, based on the findings above.
	if( nimin ) mBoundaryIMINContainer.reserve( nimin );
	if( nimax ) mBoundaryIMAXContainer.reserve( nimax );
	if( njmin ) mBoundaryJMINContainer.reserve( njmin );
	if( njmax ) mBoundaryJMAXContainer.reserve( njmax );

	// Next, loop over the markers in this zone (again) to initialize 
	// exactly as much boundary containers as we require.
	for( auto& marker: zone->GetMarker() )
	{
		// In case this marker is an interface, skip it.
		if( marker->GetTypeBC() == ETypeBC::INTERFACE ) continue;

		// Loop over each element on this marker and instantiate its boundary.
		for( auto& [index, face]: marker->GetElementFaces() )
		{
			// Based on the face type, initialzie the boundary surface.
			switch( face )
			{
				// Check the IMIN boundary.
				case(EFaceElement::IMIN): 
				{
					// Otherwise, allocate the proper boundary for this face.
					mBoundaryIMINContainer.emplace_back
					(
					 CGenericFactory::CreateBoundaryContainer(config_container, geometry_container, marker.get(), index )
					);
 
					break;
				}

				// Check the IMAX boundary.
				case(EFaceElement::IMAX): 
				{
					// Otherwise, allocate the proper boundary for this face.
					mBoundaryIMAXContainer.emplace_back
					(
					 CGenericFactory::CreateBoundaryContainer(config_container, geometry_container, marker.get(), index )
					);
 
					break;
				}

				// Check the JMIN boundary.
				case(EFaceElement::JMIN): 
				{
					// Otherwise, allocate the proper boundary for this face.
					mBoundaryJMINContainer.emplace_back
					(
					 CGenericFactory::CreateBoundaryContainer(config_container, geometry_container, marker.get(), index )
					);
 
					break;
				}

				// Check the JMAX boundary.
				case(EFaceElement::JMAX): 
				{
					// Otherwise, allocate the proper boundary for this face.
					mBoundaryJMAXContainer.emplace_back
					(
					 CGenericFactory::CreateBoundaryContainer(config_container, geometry_container, marker.get(), index )
					);
 
					break;
				}

				// Something went wrong, issue an error.
				default: ERROR("Unknown face detected.");
			}
		}
	}


	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void CEESolver::ComputeVolumeResidual
(
 as3double                  localtime,
 as3vector1d<as3double>    &monitordata,
 CPoolMatrixAS3<as3double> &workarray
)
 /*
	* Function that computes the volume terms of an EE-type PDE. 
	*/
{
	// Extract the number of integration points in 2D.
	size_t nInt2D = mStandardElementContainer->GetnInt2D();
	// Extract the quadrature integration weights on the standard element.
	auto&  wInt2D = mStandardElementContainer->GetwInt2D();

	// Borrow memory for the solution and its gradient.
	CWorkMatrixAS3<as3double> var    = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
	CWorkMatrixAS3<as3double> dVarDx = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
	CWorkMatrixAS3<as3double> dVarDy = workarray.GetWorkMatrixAS3(mNVar, nInt2D);


	// Loop over each element and compute its residual. Note, this also resets the residual.
	for( auto& physical_element: mPhysicalElementContainer )
	{
		// Extract the metrics at the volume solution points.
		auto& jac = physical_element->mMetricInt2D;
		// Reference to the current element solution.
		auto& sol = physical_element->mSol2D;
		// Reference to the current element residual.
		auto& res = physical_element->mRes2D;

		// Compute the solution and its (parametric) gradient at the volume integration points.
		mTensorProductContainer->Volume(mNVar, sol.data(),
				                            var.data(), dVarDx.data(), dVarDy.data());

		// Convert the gradient from parametric to Cartesian coordinates.
		physical_element->ConvertGradParamToCartVolInt(dVarDx, dVarDy);


		// Loop over the integration points and compute the weighted terms.
		for(size_t l=0; l<nInt2D; l++)
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
		mTensorProductContainer->ResidualVolume(mNVar, 
				                                    nullptr, 
																						dVarDx.data(), 
																						dVarDy.data(), 
																						res.data());
	}
}

//-----------------------------------------------------------------------------------

void CEESolver::ComputeSurfaceResidualIDir
(
 CGeometry                 *geometry_container,
 as3vector1d<as3double>    &monitordata,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
	* Function that computes the residual terms in the i-direction of an EE-type PDE. 
	*/
{
	// Extract pointer to the relevant grid zone.	
	auto* grid = geometry_container->GetZoneGeometry(mZoneID);

	// Extract the number of elements in this zone.
	const size_t nxElem = grid->GetnxElem();	
	const size_t nyElem = grid->GetnyElem();

	// Extract the number of integration points in 1D.
	size_t nInt1D = mStandardElementContainer->GetnInt1D();
	// Extract the quadrature integration weights in 1D on the standard element.
	auto&  wInt1D = mStandardElementContainer->GetwInt1D();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varL = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> varR = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, nInt1D);


	// First, loop over the internal faces in the i-direction and compute their residuals.
	for(size_t jElem=0; jElem<nyElem; jElem++)
	{
		for(size_t iElem=1; iElem<nxElem; iElem++)
		{
			// Deduce the indices of the left and right elements.
			const size_t IR = jElem*nxElem + iElem;
			const size_t IL = IR-1;
			
			// Get a pointer to the respective left and right element, w.r.t. this face.
			auto& elemL = mPhysicalElementContainer[IL];
			auto& elemR = mPhysicalElementContainer[IR];

			// Extract the metrics at the integration points on the imax face of the left element.
			auto& metL = elemL->mMetricIntIMax1D;
			// Reference to the left element solution.
			auto& solL = elemL->mSol2D;
			// Reference to the left element residual.
			auto& resL = elemL->mRes2D;

			// Reference to the right element solution.
			auto& solR = elemR->mSol2D;
			// Reference to the right element residual.
			auto& resR = elemR->mRes2D;

			// Compute the solution on the integration nodes of the left element.
			mTensorProductContainer->SurfaceIMAX(mNVar, solL.data(),
				                                   varL.data(), nullptr, nullptr); 

			// Compute the solution on the integration nodes of the right element.
			mTensorProductContainer->SurfaceIMIN(mNVar, solR.data(),
				                                   varR.data(), nullptr, nullptr); 

			// Compute the flux state, weighted by the integration nodes and metrics. 
			// Notice, this is based on the left state, which is the outward-pointing
			// normal vector.
			mRiemannSolverContainer->ComputeFlux(wInt1D, metL, varL, varR, flux);

			// Compute the residual on the left  element, which is on the IMAX boundary.
			mTensorProductContainer->ResidualSurfaceIMAX(mNVar, flux.data(), nullptr, nullptr, resL.data());

			// For local conservation, negate the flux, since it leaves the left element to enter the right element.
			for(size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

			// Compute the residual on the right element, which is on the IMIN boundary.
			mTensorProductContainer->ResidualSurfaceIMIN(mNVar, flux.data(), nullptr, nullptr, resR.data());
		}
	}


	// Borrow memory for the gradient of the solution, possible used in the boundary.
	CWorkMatrixAS3<as3double> dVarDxL = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> dVarDyL = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> dVarDxR = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> dVarDyR = workarray.GetWorkMatrixAS3(mNVar, nInt1D);


	// Loop over the boundary elements on the IMIN surface.
	for( auto& imin: mBoundaryIMINContainer )
	{
		// Extract the element index of this boundary face (which is on the right).
		const size_t IR = imin->GetIndex();

		// Get a pointer to the right element.
		auto& elemR = mPhysicalElementContainer[IR];

		// Extract the metrics at the integration points on the imin face of the right element.
		auto& metR = elemR->mMetricIntIMin1D;
		// Reference to the right element solution and residual.
		auto& solR = elemR->mSol2D;
		auto& resR = elemR->mRes2D;

		// Compute the solution on the integration nodes of the right element.
		mTensorProductContainer->SurfaceIMIN(mNVar, solR.data(),
			                                   varR.data(), dVarDxR.data(), dVarDyR.data());


		// Convert the gradient from parametric to Cartesian coordinates.
		elemR->ConvertGradParamToCartSurfIMinInt(dVarDxR, dVarDyR);

		// Compute the boundary state.
		imin->ComputeBoundaryState(metR, 
				                       varR, dVarDxR, dVarDyR, 
															 varL, dVarDxL, dVarDyL);

		// Compute the flux state, weighted by the integration nodes and metrics. 
		// Notice, this is taken from the perspective of the right state, which 
		// is the inward-pointing normal vector.
		mRiemannSolverContainer->ComputeFlux(wInt1D, metR, varL, varR, flux);

		// Compute the residual on the right element, which is on the IMIN boundary.
		mTensorProductContainer->ResidualSurfaceIMIN(mNVar, flux.data(), nullptr, nullptr, resR.data());
	}


	// Loop over the boundary elements on the IMAX surface.
	for( auto& imax: mBoundaryIMAXContainer )
	{
		// Extract the element index of this boundary face (which is on the left).
		const size_t IL = imax->GetIndex();

		// Get a pointer to the left element.
		auto& elemL = mPhysicalElementContainer[IL];

		// Extract the metrics at the integration points on the imax face of the left element.
		auto& metL = elemL->mMetricIntIMax1D;
		// Reference to the left element solution and residual.
		auto& solL = elemL->mSol2D;
		auto& resL = elemL->mRes2D;

		// Compute the solution on the integration nodes of the left element.
		mTensorProductContainer->SurfaceIMIN(mNVar, solL.data(),
			                                   varL.data(), dVarDxL.data(), dVarDyL.data()); 

		// Convert the gradient from parametric to Cartesian coordinates.
		elemL->ConvertGradParamToCartSurfIMaxInt(dVarDxL, dVarDyL);

		// Compute the boundary state.
		imax->ComputeBoundaryState(metL, 
				                       varL, dVarDxL, dVarDyL,
															 varR, dVarDxR, dVarDyR);

		// Compute the flux state, weighted by the integration nodes and metrics. 
		// Notice, this is taken from the perspective of the left state, which 
		// is the outward-pointing normal vector.
		mRiemannSolverContainer->ComputeFlux(wInt1D, metL, varL, varR, flux);

		// Compute the residual on the left element, which is on the IMAX boundary.
		mTensorProductContainer->ResidualSurfaceIMAX(mNVar, flux.data(), nullptr, nullptr, resL.data());
	}
}

//-----------------------------------------------------------------------------------

void CEESolver::ComputeSurfaceResidualJDir
(
 CGeometry                 *geometry_container,
 as3vector1d<as3double>    &monitordata,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
	* Function that computes the residual terms in the j-direction of an EE-type PDE. 
	*/
{
	// Extract pointer to the relevant grid zone.	
	auto* grid = geometry_container->GetZoneGeometry(mZoneID);

	// Extract the number of elements in this zone.
	const size_t nxElem = grid->GetnxElem();	
	const size_t nyElem = grid->GetnyElem();

	// Extract the number of integration points in 1D.
	size_t nInt1D = mStandardElementContainer->GetnInt1D();
	// Extract the quadrature integration weights in 1D on the standard element.
	auto&  wInt1D = mStandardElementContainer->GetwInt1D();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varB = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> varT = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, nInt1D);


	// First, loop over the internal faces in the j-direction and compute their residuals.
	for(size_t jElem=1; jElem<nyElem; jElem++)
	{
		for(size_t iElem=0; iElem<nxElem; iElem++)
		{
			// Deduce the indices of the left and right elements.
			const size_t IT = jElem*nxElem + iElem;
			const size_t IB = IT-nxElem;
			
			// Get a pointer to the respective bottom and top element, w.r.t. this face.
			auto& elemB = mPhysicalElementContainer[IB];
			auto& elemT = mPhysicalElementContainer[IT];

			// Extract the metrics at the integration points on the jmax face of the bottom element.
			auto& metB = elemB->mMetricIntJMax1D;
			// Reference to the bottom element solution.
			auto& solB = elemB->mSol2D;
			// Reference to the bottom element residual.
			auto& resB = elemB->mRes2D;

			// Reference to the top element solution.
			auto& solT = elemT->mSol2D;
			// Reference to the top element residual.
			auto& resT = elemT->mRes2D;

			// Compute the solution on the integration nodes of the bottom element.
			mTensorProductContainer->SurfaceJMAX(mNVar, solB.data(),
				                                   varB.data(), nullptr, nullptr); 

			// Compute the solution on the integration nodes of the top element.
			mTensorProductContainer->SurfaceJMIN(mNVar, solT.data(),
				                                   varT.data(), nullptr, nullptr); 

			// Compute the flux state, weighted by the integration nodes and metrics. 
			// Notice, this is based on the bottom state, which is the outward-pointing
			// normal vector.
			mRiemannSolverContainer->ComputeFlux(wInt1D, metB, varB, varT, flux);

			// Compute the residual on the bottom  element, which is on the JMAX boundary.
			mTensorProductContainer->ResidualSurfaceJMAX(mNVar, flux.data(), nullptr, nullptr, resB.data());

			// For local conservation, negate the flux, since it leaves the bottom element to enter the top element.
			for(size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

			// Compute the residual on the top element, which is on the JMIN boundary.
			mTensorProductContainer->ResidualSurfaceJMIN(mNVar, flux.data(), nullptr, nullptr, resT.data());
		}
	}


	// Borrow memory for the gradient of the solution, possible used in the boundary.
	CWorkMatrixAS3<as3double> dVarDxB = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> dVarDyB = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> dVarDxT = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> dVarDyT = workarray.GetWorkMatrixAS3(mNVar, nInt1D);


	// Loop over the boundary elements on the JMIN surface.
	for( auto& jmin: mBoundaryJMINContainer )
	{
		// Extract the element index of this boundary face (which is on the top).
		const size_t IT = jmin->GetIndex();

		// Get a pointer to the top element.
		auto& elemT = mPhysicalElementContainer[IT];

		// Extract the metrics at the integration points on the jmin face of the top element.
		auto& metT = elemT->mMetricIntJMin1D;
		// Reference to the right element solution and residual.
		auto& solT = elemT->mSol2D;
		auto& resT = elemT->mRes2D;

		// Compute the solution on the integration nodes of the top element.
		mTensorProductContainer->SurfaceJMIN(mNVar, solT.data(),
			                                   varT.data(), dVarDxT.data(), dVarDyT.data());


		// Convert the gradient from parametric to Cartesian coordinates.
		elemT->ConvertGradParamToCartSurfJMinInt(dVarDxT, dVarDyT);

		// Compute the boundary state.
		jmin->ComputeBoundaryState(metT, 
				                       varT, dVarDxT, dVarDyT, 
															 varB, dVarDxB, dVarDyB);

		// Compute the flux state, weighted by the integration nodes and metrics. 
		// Notice, this is taken from the perspective of the top state, which 
		// is the inward-pointing normal vector.
		mRiemannSolverContainer->ComputeFlux(wInt1D, metT, varB, varT, flux);

		// Compute the residual on the top element, which is on the JMIN boundary.
		mTensorProductContainer->ResidualSurfaceJMIN(mNVar, flux.data(), nullptr, nullptr, resT.data());
	}


	// Loop over the boundary elements on the JMAX surface.
	for( auto& jmax: mBoundaryJMAXContainer )
	{
		// Extract the element index of this boundary face (which is on the bottom).
		const size_t IB = jmax->GetIndex();

		// Get a pointer to the bottom element.
		auto& elemB = mPhysicalElementContainer[IB];

		// Extract the metrics at the integration points on the jmax face of the bottom element.
		auto& metB = elemB->mMetricIntJMax1D;
		// Reference to the bottom element solution and residual.
		auto& solB = elemB->mSol2D;
		auto& resB = elemB->mRes2D;

		// Compute the solution on the integration nodes of the bottom element.
		mTensorProductContainer->SurfaceJMIN(mNVar, solB.data(),
			                                   varB.data(), dVarDxB.data(), dVarDyB.data()); 

		// Convert the gradient from parametric to Cartesian coordinates.
		elemB->ConvertGradParamToCartSurfJMaxInt(dVarDxB, dVarDyB);

		// Compute the boundary state.
		jmax->ComputeBoundaryState(metB, 
				                       varB, dVarDxB, dVarDyB,
															 varT, dVarDxT, dVarDyT);

		// Compute the flux state, weighted by the integration nodes and metrics. 
		// Notice, this is taken from the perspective of the bottom state, which 
		// is the outward-pointing normal vector.
		mRiemannSolverContainer->ComputeFlux(wInt1D, metB, varB, varT, flux);

		// Compute the residual on the bottom element, which is on the JMAX boundary.
		mTensorProductContainer->ResidualSurfaceJMAX(mNVar, flux.data(), nullptr, nullptr, resB.data());
	}
}






