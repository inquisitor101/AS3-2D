#include "iteration_structure.hpp"


//-----------------------------------------------------------------------------------
// CIteration member functions.
//-----------------------------------------------------------------------------------


CIteration::CIteration
(
 CConfig                               *config_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
 /*
	* Constructor for the iteration class.
	*/
{
	// Check the relevant number of data entries required in the work array.
	for(unsigned short iZone=0; iZone<config_container->GetnZone(); iZone++)
	{
		// For now, take the maximum number of items.
		switch( config_container->GetTypeSolver(iZone) )
		{
			case(ETypeSolver::EE):
			{
				const size_t nItem2D = 3;  // volume  terms needed.
				const size_t nItem1D = 2;  // surface terms needed.
				
				const size_t nVar    = solver_container[iZone]->GetnVar();
				const size_t nInt1D  = solver_container[iZone]->GetStandardElement()->GetnInt1D();
				const size_t nInt2D  = solver_container[iZone]->GetStandardElement()->GetnInt2D();

				// Compute the total number of required volume and surface terms in the work array.
				const size_t nVol  = nItem2D*nInt2D*nVar;
				const size_t nSurf = nItem1D*nInt1D*nVar;

				// Take the maximum storage between the volume and surface terms.
				const size_t nData = std::max( nVol, nSurf );

				// Take whichever is the max possible storage across all zones (can be inefficient).
				mNWorkData = std::max( mNWorkData, nData );
				
				break;
			}

			default: ERROR("Unknown solver type.");
		}
	}
}

//-----------------------------------------------------------------------------------

CIteration::~CIteration
(
 void
)
 /*
	* Destructor, which cleans up after the iteration class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::PreProcessIteration
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3double                                 localtime 
)
 /*
	* Function that preprocesses the solution, before sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::PostProcessIteration
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3double                                 localtime 
)
 /*
	* Function that postprocesses the solution, after sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::ComputeResidual
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3double                                 localtime
)
 /*
	* Function that computes the residual in all zones.
	*/
{
	// Get the total number of elements in all zones.
	const size_t nElemTotal = openmp_container->GetnIndexVolume();
	// Get the total number of internal i-faces in all zones.
	const size_t nInternalIFace = openmp_container->GetnInternIFace();
	// Get the total number of internal j-faces in all zones.
	const size_t nInternalJFace = openmp_container->GetnInternJFace();

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
		// Deduce the current element's zone and index.
		const unsigned short iZone = openmp_container->GetIndexVolume(i)->mZone;
		const unsigned int   iElem = openmp_container->GetIndexVolume(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the relevant grid.
		auto* grid    = geometry_container->GetZoneGeometry(iZone);

		// Compute the volume terms on this element. Note, this step also initializes the residual.
		solver->ComputeVolumeResidual(grid, workarray, localtime, iElem); 
	}


	// TESTING: start
	//----------------------------------------------------------------
	const auto& idir_faces_multizone = geometry_container->GetMultizoneIFaces();
	const auto& jdir_faces_multizone = geometry_container->GetMultizoneJFaces();
	
  const size_t nFacesJDir = jdir_faces_multizone.GetnFacesTotal(); 
  const size_t nFacesIDir = idir_faces_multizone.GetnFacesTotal(); 
	
  // NOTE, mResMinus is always known how much to allocate in each zone, since it is the convention
	// that, regardless if the face is internal, boundary or interface, mResMin are always book-keeped
	// to avoid race conditions. Hence, can allocate them per zone without additional information.

	// TODO: move the compute residuals to a separate numerics_container and pass each solver to it.

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesIDir; i++)
	{
		// Get current face type.
		const auto face_type = idir_faces_multizone.GetFaceTypeFromIndex(i);

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the local index.
        const size_t iFaceLocal = idir_faces_multizone.GetIndexInternalFace(i);
			
        // Extract the face information, essentially its indices.
        const auto face_info = idir_faces_multizone.GetInternalElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
       
        // Extract the current face.
				const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );
				
        const unsigned short iZone  = family_face.GetiZone();
				const size_t         iElemL = internal_face.mIndexElementM;

				const auto* grid = geometry_container->GetZoneGeometry(iZone);

				solver_container[iZone]->ComputeSurfaceResidualIDir_NEW(grid, workarray, localtime, iElemL);

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the local index.
        const size_t iFaceLocal = idir_faces_multizone.GetIndexInterfaceFace(i);
			
        // Extract the face information, essentially its indices.
        const auto face_info = idir_faces_multizone.GetInterfaceElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInterfaceFacesFamily( face_info.mIndexFamily );
       
        // TESTING: doesn't matter the index 0 or not, as it doesnt use any member variables of IInterface.
        interface_container[0]->ComputeInterfaceResidual_NEW(solver_container, family_face, face_info, workarray, localtime);
        break;
      }

			default: ERROR("Unknown face type, could not compute its residual.");
		}

	}

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesIDir; i++)
	{
		// Get current face type.
		const auto face_type = idir_faces_multizone.GetFaceTypeFromIndex(i);

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the local index.
        const size_t iFaceLocal = idir_faces_multizone.GetIndexInternalFace(i);

        // Extract the face information, essentially its indices.
        const auto face_info = idir_faces_multizone.GetInternalElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
        
        // Extract the current face.
        const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );

				const unsigned short iZone  = family_face.GetiZone();
				const size_t         iElemL = internal_face.mIndexElementM;

        auto* physical_element = solver_container[iZone]->GetPhysicalElement(iElemL);

        auto& tmpL = physical_element->mResMinus;
        auto& resL = physical_element->mRes2D;
        
        for(size_t l=0; l<resL.size(); l++) resL[l] += tmpL[l];

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the local index.
        const size_t iFaceLocal = idir_faces_multizone.GetIndexInterfaceFace(i);

        // Extract the face information, essentially its indices.
        const auto face_info = idir_faces_multizone.GetInterfaceElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInterfaceFacesFamily( face_info.mIndexFamily );

        // Extract the actual face.
        const auto& interface_face = family_face.GetInterfaceFace( face_info.mIndexFace ); 

        const unsigned short iZoneM = family_face.GetiZone();
        const unsigned short iZoneP = family_face.GetjZone();

        const size_t iElemM = interface_face.mIndexElementI;
        const size_t iElemP = interface_face.mIndexElementJ;

        const EFaceLocation iFaceM = family_face.GetiFaceLocation();
        const EFaceLocation iFaceP = family_face.GetjFaceLocation();

        auto* physical_element_m = solver_container[iZoneM]->GetPhysicalElement(iElemM);
        auto* physical_element_p = solver_container[iZoneP]->GetPhysicalElement(iElemP);

        auto& tmp_m = physical_element_m->mResMinus;
        auto& res_m = physical_element_m->mRes2D;

        auto& tmp_p = physical_element_p->mResMinus;
        auto& res_p = physical_element_p->mRes2D;

        if( iFaceM == EFaceLocation::IMAX || iFaceM == EFaceLocation::JMAX )
        {
          for(size_t l=0; l<res_m.size(); l++) res_m[l] += tmp_m[l];
        }

        if( iFaceP == EFaceLocation::IMAX || iFaceP == EFaceLocation::JMAX )
        {
          for(size_t l=0; l<res_p.size(); l++) res_p[l] += tmp_p[l];
        }

        break;
      }

			default: ERROR("Unknown face type, could not compute its residual.");
		}
	}






#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesJDir; i++)
	{
		// Get current face type.
		const auto face_type = jdir_faces_multizone.GetFaceTypeFromIndex(i);

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the local index.
        const size_t iFaceLocal = jdir_faces_multizone.GetIndexInternalFace(i);
				
        // Extract the face information, essentially its indices.
        const auto face_info = jdir_faces_multizone.GetInternalElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
       
        // Extract the current face.
        const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );

        const unsigned short iZone  = family_face.GetiZone();
        const size_t         iElemB = internal_face.mIndexElementM;

        const auto* grid = geometry_container->GetZoneGeometry(iZone);

				solver_container[iZone]->ComputeSurfaceResidualJDir_NEW(grid, workarray, localtime, iElemB);

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the local index.
        const size_t iFaceLocal = jdir_faces_multizone.GetIndexInterfaceFace(i);
        
        // Extract the face information, essentially its indices.
        const auto face_info = jdir_faces_multizone.GetInterfaceElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInterfaceFacesFamily( face_info.mIndexFamily );

        // TESTING: doesn't matter the index 0 or not, as it doesnt use any member variables of IInterface.
        interface_container[0]->ComputeInterfaceResidual_NEW(solver_container, family_face, face_info, workarray, localtime);
        break;
      }

			default: ERROR("Unknown face type, could not compute its residual.");
		}
	}

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesJDir; i++)
	{
		// Get current face type.
		const auto face_type = jdir_faces_multizone.GetFaceTypeFromIndex(i);

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the local index.
        const size_t iFaceLocal = jdir_faces_multizone.GetIndexInternalFace(i);
				
        // Extract the face information, essentially its indices.
        const auto face_info = jdir_faces_multizone.GetInternalElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
       
        // Extract the current face.
        const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );

				const unsigned short iZone  = family_face.GetiZone();
				const size_t         iElemB = internal_face.mIndexElementM;

        auto* physical_element = solver_container[iZone]->GetPhysicalElement(iElemB);

        auto& tmpB = physical_element->mResMinus;
        auto& resB = physical_element->mRes2D;
        
        for(size_t l=0; l<resB.size(); l++) resB[l] += tmpB[l];

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the local index.
        const size_t iFaceLocal = jdir_faces_multizone.GetIndexInterfaceFace(i);
        
        // Extract the face information, essentially its indices.
        const auto face_info = jdir_faces_multizone.GetInterfaceElementFaceIndex(iFaceLocal);

        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInterfaceFacesFamily( face_info.mIndexFamily );

        // Extract the actual face.
        const auto& interface_face = family_face.GetInterfaceFace( face_info.mIndexFace ); 

        const unsigned short iZoneM = family_face.GetiZone();
        const unsigned short iZoneP = family_face.GetjZone();

        const size_t iElemM = interface_face.mIndexElementI;
        const size_t iElemP = interface_face.mIndexElementJ;

        const EFaceLocation iFaceM = family_face.GetiFaceLocation();
        const EFaceLocation iFaceP = family_face.GetjFaceLocation();

        auto* physical_element_m = solver_container[iZoneM]->GetPhysicalElement(iElemM);
        auto* physical_element_p = solver_container[iZoneP]->GetPhysicalElement(iElemP);

        auto& tmp_m = physical_element_m->mResMinus;
        auto& res_m = physical_element_m->mRes2D;

        auto& tmp_p = physical_element_p->mResMinus;
        auto& res_p = physical_element_p->mRes2D;

        if( iFaceM == EFaceLocation::IMAX || iFaceM == EFaceLocation::JMAX )
        {
          for(size_t l=0; l<res_m.size(); l++) res_m[l] += tmp_m[l];
        }

        if( iFaceP == EFaceLocation::IMAX || iFaceP == EFaceLocation::JMAX )
        {
          for(size_t l=0; l<res_p.size(); l++) res_p[l] += tmp_p[l];
        }

        break;
      }

			default: ERROR("Unknown face type, could not compute its residual.");
		}
	}







	//// Loop over the interface boundaries, if need be.
	//for( auto& interface: interface_container )
	//{
	//	interface->ComputeInterfaceResidual(solver_container, workarray);
	//}
	//----------------------------------------------------------------
	// TESTING: over




	// Multiply by the inverse mass matrix.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
		// Deduce the current element's zone and index.
		const unsigned short iZone = openmp_container->GetIndexVolume(i)->mZone;
		const unsigned int   iElem = openmp_container->GetIndexVolume(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = solver_container[iZone];
		// Extract the relevant physical element.
		auto* element = solver->GetPhysicalElement(iElem);

		// Get the inverse of the mass matrix.
		auto& m = element->mInvMassMatrix;
		// Get the residual on this element.
		auto& r = element->mRes2D;

		// Perform a matrix-matrix multiplication to obtain the residual.
		NLinearAlgebra::MatrixVectorTransMult(m, r, r);
	}
}

//-----------------------------------------------------------------------------------

void CIteration::GridSweep
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 COpenMP                                  *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container, 
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 as3double                                 localtime 
)
 /*
	* Function that performs a grid sweep over all the zones. 
	*/
{
	// Initialize a work array, to avoid multiple memory allocations.
	// Note, during parallelization, this needs to be allocated inside 
	// the (shared memory) parallel region -- not here.
	CPoolMatrixAS3<as3double> workarray(mNWorkData);


	// Check for any preprocessing steps.
	PreProcessIteration(config_container,
			                geometry_container,
											openmp_container,
											solver_container,
											interface_container,
											workarray,
											localtime);


	// Compute the residual over all zones.
	ComputeResidual(config_container,
			            geometry_container,
									openmp_container,
									solver_container,
									interface_container,
									workarray,
									localtime);


	// Check for any postprocessing steps.
	PostProcessIteration(config_container,
			                 geometry_container,
											 openmp_container,
											 solver_container,
											 interface_container,
											 workarray,
											 localtime);
}


