#include "interface_structure.hpp"


//-----------------------------------------------------------------------------------
// IInterface member functions.
//-----------------------------------------------------------------------------------


IInterface::IInterface
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 CInterfaceParamMarker                 *param_container,
 const CMarker                         *imarker_container,
 const CMarker                         *jmarker_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
	:
		mIName( imarker_container->GetNameMarker() ),
		mJName( jmarker_container->GetNameMarker() ),
		mIZone( imarker_container->GetZoneID()     ),
		mJZone( jmarker_container->GetZoneID()     )
 /*
	* Constructor for the interface zone interface class.
	*/
{
	// Temporary lambda to find the face direction on a given marker.
	auto lFaceDir = [=](const CMarker *marker_container) -> EFaceElement
	{
		// Select the face direction based on the first element index.
		EFaceElement iface = marker_container->GetElementFaces(0).mFace; 
		
		// Loop over each element and ensure the face is constant.
		for( auto& marker: marker_container->GetElementFaces() )
		{
			if( marker.mFace != iface )
			{
				ERROR(marker_container->GetNameMarker() + " must have the same face direction.");
			}
		}

		// Return the face direction.
		return iface;
	};

	// Deduce the faces of each marker. Note, it suffices to consider only the 
	// first element, as the entire marker must have the same face direction.
	mIFace = lFaceDir(imarker_container); 
	mJFace = lFaceDir(jmarker_container); 

	// Deduce the number of elements on both markers.
	mNElem = imarker_container->GetnElem();
	
	// Check that the number of elements is not zero.
	if( mNElem == 0 ) ERROR("Interface markers must not be empty.");
	
	// Ensure that both markers have the same number of elements. 
	if( mNElem != jmarker_container->GetnElem() )
	{
		ERROR("Interface markers must share the same number of elements.");
	}

	// Initialize the elements of both pair of markers.
	mIndexElement.reserve(mNElem);

	for(unsigned int i=0; i<mNElem; i++)
	{
		// The assumption in AS3 is that the matching face is reversed. 
		// This is because all zones use a clockwise convention to tag 
		// their boundary markers. The the matching (j-)index is:
		const unsigned int j = mNElem - i - 1;

		// Deduce the actual element indices, not their marker index.
		const unsigned int I = imarker_container->GetElementFaces(i).mIndex;
		const unsigned int J = jmarker_container->GetElementFaces(j).mIndex;

		mIndexElement.emplace_back(I,J);
	}

	// Process the markers, to ensure they coincide geometrically.
	ProcessMatchingMarkers(config_container, 
			                   geometry_container, 
												 imarker_container, 
												 jmarker_container,
												 param_container);
}

//-----------------------------------------------------------------------------------

IInterface::~IInterface
(
 void
)
 /*
	* Destructor, which cleans up after the interface zone interface class.
	*/
{

}

//-----------------------------------------------------------------------------------

void IInterface::ProcessMatchingMarkers
(
 CConfig               *config_container,
 CGeometry             *geometry_container,
 const CMarker         *imarker_container,
 const CMarker         *jmarker_container,
 CInterfaceParamMarker *param_container
)
 /*
	* Function that processes each pair of markers, such that their common face coincides.
	*/
{
	// Extract the grid geometry in each of these zones.
	auto* igrid = geometry_container->GetZoneGeometry(mIZone);
	auto* jgrid = geometry_container->GetZoneGeometry(mJZone);

	// Extract the properties of each marker region (element indices and faces).
	auto& imarker = imarker_container->GetElementFaces(); 
	auto& jmarker = jmarker_container->GetElementFaces();

	// Ensure the number of indices in each marker matches.
	if( (imarker.size() != jmarker.size()) || (imarker.size() != mNElem) ) 
	{
		ERROR("Interface markers have different number of elements.");
	}

	// Loop over each pair of elements and check that their faces coincide.
	for(size_t i=0; i<mNElem; i++)
	{
		// The assumption in AS3 is that the matching face is reversed. 
		// This is because all zones use a clockwise convention to tag 
		// their boundary markers. The the matching (j-)index is:
		const size_t j = mNElem - i - 1;

		// Extract the nodal indices of the first element on this marker.
		auto& icoor = igrid->GetElementGeometry( imarker[i].mIndex )->GetCoordSolDOFs(); 
		auto& jcoor = jgrid->GetElementGeometry( jmarker[j].mIndex )->GetCoordSolDOFs(); 
		
		// Extract the nodal indices from the face type.
		auto& inode = igrid->GetFaceNodalIndices( imarker[i].mFace );
		auto& jnode = jgrid->GetFaceNodalIndices( jmarker[j].mFace ); 

		// NOTE 
		// For now, force the polynomial orders to be the same. Otherwise, we 
		// have to come up with an integration rule that is based on the higher
		// order polynomial, to ensure local conservation.
		if( inode.size() != jnode.size() )
		{
			ERROR("Polynomial order must be identical (for now) in zones: "
					  + std::to_string(mIZone) + ", " + std::to_string(mJZone));
		}

		// Relative tolerance value.
		const as3double tol = static_cast<as3double>( 1.0e-8 );
		
		// Check if the periodic faces coincide, after translation.
		for(size_t k=0; k<inode.size(); k++)
		{
			// Extract coordinates for the current  marker.
			const as3double ix = icoor(0, inode[k]);
			const as3double iy = icoor(1, inode[k]);
			// Extract coordinates for the matching marker.
			const as3double jx = jcoor(0, jnode[k]);
			const as3double jy = jcoor(1, jnode[k]);
			
			// Compute the difference between them in absolute.
			const as3double dx = std::abs( ix - jx + param_container->mVectorTrans[0] );
			const as3double dy = std::abs( iy - jy + param_container->mVectorTrans[1] );

			// Relative error, based on the values of the coordinates.
			const as3double xtol = tol*std::max( std::abs(ix), std::abs(jx) );
			const as3double ytol = tol*std::max( std::abs(iy), std::abs(jy) );

			// If the boundaries do not coincide, abort with an error.
			if( (dx > xtol) || (dy > ytol) ) 
			{
				ERROR("Interface boundaries: " + mIName + ", " + mJName + " do not match.");
			}
		}
	}
}


//-----------------------------------------------------------------------------------
// CEEInterface member functions.
//-----------------------------------------------------------------------------------


CEEInterface::CEEInterface
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 CInterfaceParamMarker                 *param_container,
 const CMarker                         *imarker_container,
 const CMarker                         *jmarker_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
	:
		IInterface(config_container, 
				       geometry_container, 
							 param_container,
							 imarker_container,
							 jmarker_container,
				       solver_container)
 /*
	* Constructor for the (non-linear) Euler equations zone interface class.
	*/
{
	// Check that the (owner) imarker indeed is of type interface.
	if( imarker_container->GetTypeBC() != ETypeBC::INTERFACE )
	{
		ERROR(imarker_container->GetNameMarker() + " must be an interface boundary.");
	}

	// Check that the (matching) jmarker indeed is of type interface.
	if( jmarker_container->GetTypeBC() != ETypeBC::INTERFACE )
	{
		ERROR(jmarker_container->GetNameMarker() + " must be an interface boundary.");
	}


	// Extract the polynomial orders in each marker zone.
	unsigned short ipoly = solver_container[mIZone]->GetStandardElement()->GetnPolySol();
	unsigned short jpoly = solver_container[mJZone]->GetStandardElement()->GetnPolySol();
	
	// Extract the number of integration points in each marker zone.
	unsigned short inint = solver_container[mIZone]->GetStandardElement()->GetnInt1D();
	unsigned short jnint = solver_container[mJZone]->GetStandardElement()->GetnInt1D();

	// Extract the type of DOFs in each zone.
	ETypeDOF idof = solver_container[mIZone]->GetStandardElement()->GetTypeDOFsSol();
	ETypeDOF jdof = solver_container[mJZone]->GetStandardElement()->GetTypeDOFsSol();

	// Only choose EQD when both markers are EQD. Otherwise, always select LGL.
	ETypeDOF dof = ((idof == ETypeDOF::EQD) && (jdof == ETypeDOF::EQD)) ? ETypeDOF::EQD : ETypeDOF::LGL;

	// Take the integration rule based on the highest polynomial.
	mNInt1D = std::max( inint, jnint );

	// Instantiate the appropriate (temporary) standard element containers in each zone.
	CStandardElement ielement(dof, ipoly, mNInt1D);
	CStandardElement jelement(dof, jpoly, mNInt1D);

	// Obtain the integration weights on this interface.
	mWInt1D = ielement.GetwInt1D();

	// Instantiate the tensor-product containers in iZone and jZone.
	mITensorProductContainer = CGenericFactory::CreateTensorContainer( &ielement );
	mJTensorProductContainer = CGenericFactory::CreateTensorContainer( &jelement );

	// Extract the type of Riemann solver in each zone.
	auto iriemann = solver_container[mIZone]->GetRiemannSolver()->GetTypeRiemannSolver();  
	auto jriemann = solver_container[mJZone]->GetRiemannSolver()->GetTypeRiemannSolver();
	
	// For now, force the Riemann solvers in both zones to be identical.
	if( iriemann != jriemann )
	{
		ERROR("For now, Riemann solvers at an interface must be identical.");
	}

	// Instantiate a Riemann solver on this interface.
	mRiemannSolverContainer = CGenericFactory::CreateRiemannSolverContainer( config_container, iriemann );

	// Ensure the number of working variables in the iZone is as expected.
	if( mNVar != solver_container[mIZone]->GetnVar() )
	{
		ERROR("Number of variables mismatches in iZone: " + std::to_string(mIZone));
	}
	if( mNVar != solver_container[mJZone]->GetnVar() )
	{
		ERROR("Number of variables mismatches in jZone: " + std::to_string(mJZone));
	}
}

//-----------------------------------------------------------------------------------

CEEInterface::~CEEInterface
(
 void
)
 /*
	* Destructor, which cleans up after the Euler equations zone interface class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEEInterface::ComputeInterfaceResidual
(
 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
 CPoolMatrixAS3<as3double>             &workarray
)
 /*
	* Function that computes the residual on the zone interface marker. 
	*/
{
	// Get the appropriate functions for interpolation.
	auto InterpSurfaceI  = GetFuncPointerInterpFaceI();
	auto InterpSurfaceJ  = GetFuncPointerInterpFaceJ();

	// Get the appropriate functions for residual contribution.
	auto ComputeResFaceI = GetFuncPointerResidualFaceI();
	auto ComputeResFaceJ = GetFuncPointerResidualFaceJ();

	// Get the solvers of this class.
	auto& isolver = solver_container[mIZone];
	auto& jsolver = solver_container[mJZone];

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varI = workarray.GetWorkMatrixAS3(mNVar, mNInt1D);
	CWorkMatrixAS3<as3double> varJ = workarray.GetWorkMatrixAS3(mNVar, mNInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, mNInt1D);


	// Loop over the owned elements on this interface.
	for( auto& [I, J]: mIndexElement )
	{
		// Get a pointer to the respective elements sharing this face.
		auto* ielem = isolver->GetPhysicalElement(I);
		auto* jelem = jsolver->GetPhysicalElement(J);

		// Extract the metrics at the integration points on the owner face.
		auto& metI = ielem->GetSurfaceMetricInt(mIFace);
		// Reference to the owner element solution.
		auto& solI = ielem->mSol2D;
		// Reference to the owner element residual.
		auto& resI = ielem->mRes2D;

		// Reference to the matching element solution.
		auto& solJ = jelem->mSol2D;
		// Reference to the matching element residual.
		auto& resJ = jelem->mRes2D;

		// Compute the solution on the integration nodes of the owner element.
		InterpSurfaceI(mNVar, solI.data(), varI.data(), nullptr, nullptr); 

		// Compute the solution on the integration nodes of the matching element.
		InterpSurfaceJ(mNVar, solJ.data(), varJ.data(), nullptr, nullptr); 

		// Compute the flux state, weighted by the integration nodes and metrics. 
		// Notice, this is based on the owner state.
		mRiemannSolverContainer->ComputeFlux(mWInt1D, metI, varI, varJ, flux);

		// Compute the residual on the owned element, which is on the iface boundary.
		ComputeResFaceI(mNVar, flux.data(), nullptr, nullptr, resI.data());

		// For local conservation, negate the flux, since it leaves the owner element 
		// to enter the matching element.
		for(size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

		// Compute the residual on the matching element, which is on the jface boundary.
		ComputeResFaceJ(mNVar, flux.data(), nullptr, nullptr, resJ.data());
	}
}








