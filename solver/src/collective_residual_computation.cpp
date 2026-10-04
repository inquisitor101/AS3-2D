#include "collective_residual_computation.hpp"


//-----------------------------------------------------------------------------------
// NResidualComputation namespace functions.
//-----------------------------------------------------------------------------------

void NResidualComputation::ComputeVolumeResidualsCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  *
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

	// Get the total number of elements in all zones.
	const auto nElemTotal = multizone_geometry_container->GetnElemTotal();

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
		// Extract the element indices.
		const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);

		// Deduce the current element's zone and index.
		const auto iZone = elem_info.mIndexZone; 
		const auto iElem = elem_info.mIndexElem; 

		// Extract the relevant solver.
		auto* solver = multizone_solver_container->GetSinglezoneSolver(iZone);
		// Extract the relevant grid.
		auto* grid   = multizone_geometry_container->GetSinglezoneGeometry(iZone);

		// Compute the volume terms on this element. Note, this step also initializes the residual.
		solver->ComputeVolumeResidual(grid, workarray, localtime, iElem); 
	}
}

//-----------------------------------------------------------------------------------

void NResidualComputation::ComputeIFaceResidualsCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  *
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

	// Get the total faces in the i-direction.
	const auto& idir_faces_multizone = multizone_geometry_container->GetMultizoneIFaces();
	
	// Get the total number faces in the i--direction.
  const auto nFacesIDir = idir_faces_multizone.GetnFacesTotal(); 

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesIDir; i++)
	{
		// Get the relevant flattened index of the current face.
    const auto face_info = multizone_geometry_container->GetFlattenedIndexIFaceLoadBalanced(i);
    // Extract the current type of face.
    const auto face_type = face_info.mFaceType;

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
        // Extract the current face.
				const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );
			
        const auto iZone  = family_face.GetiZone();
				const auto iElemL = internal_face.mIndexElementM;

				const auto* grid = multizone_geometry_container->GetSinglezoneGeometry(iZone);

				// Compute the internal face residual.
				multizone_solver_container->GetSinglezoneSolver(iZone)->ComputeSurfaceResidualIDir(grid, workarray, localtime, iElemL);

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the relevant family index, which is the same as the interface container's index.
        const auto iFamily = face_info.mIndexFamily;
        // Extract the associated interface index.
        const auto iInterface = multizone_solver_container->GetIndexInterfaceIFace(iFamily); 
        
        // Extract the associated interface object.
        const auto* interface = multizone_solver_container->GetSinglezoneInterface(iInterface); 
        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInterfaceFacesFamily( iFamily );
       
        // Compute the interface residual.
        interface->ComputeInterfaceResidual(multizone_solver_container, family_face, face_info, workarray, localtime);
        
        break;
      }

			default: ERROR("Unknown IFace type, could not compute its residual.");
		}
	}
}

//-----------------------------------------------------------------------------------

void NResidualComputation::ComputeJFaceResidualsCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  *
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

	// Get the total faces in the j-direction.
	const auto& jdir_faces_multizone = multizone_geometry_container->GetMultizoneJFaces();
	
	// Get the total number faces in the  j-direction.
  const auto nFacesJDir = jdir_faces_multizone.GetnFacesTotal(); 

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesJDir; i++)
	{
		// Get the relevant flattened index of the current face.
    //const auto face_info = multizone_geometry_container->GetFlattenedIndexJFace(i);
    const auto face_info = multizone_geometry_container->GetFlattenedIndexJFaceLoadBalanced(i);
    // Extract the current type of face.
    const auto face_type = face_info.mFaceType;

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
       
        // Extract the current face.
        const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );

        const auto iZone  = family_face.GetiZone();
        const auto iElemB = internal_face.mIndexElementM;

        const auto* grid = multizone_geometry_container->GetSinglezoneGeometry(iZone);

				// Compute the internal face residual. 
				multizone_solver_container->GetSinglezoneSolver(iZone)->ComputeSurfaceResidualJDir(grid, workarray, localtime, iElemB);

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the relevant family index, which is the same as the interface container's index.
        const auto iFamily = face_info.mIndexFamily;
        // Extract the associated interface index.
        const auto iInterface = multizone_solver_container->GetIndexInterfaceJFace(iFamily); 
        
        // Extract the associated interface object.
        const auto* interface = multizone_solver_container->GetSinglezoneInterface(iInterface); 
        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInterfaceFacesFamily( iFamily );

        // Compute the interface residual.
        interface->ComputeInterfaceResidual(multizone_solver_container, family_face, face_info, workarray, localtime);
        
        break;
      }

			default: ERROR("Unknown JFace type, could not compute its residual.");
		}
	}
}

//-----------------------------------------------------------------------------------

void NResidualComputation::AccumulateIFaceResidualsCollective
(
 const CMultizoneGeometry *multizone_geometry_container,
 CMultizoneSolver         *multizone_solver_container
)
 /*
  *
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

	// Get the total faces in the i-direction.
	const auto& idir_faces_multizone = multizone_geometry_container->GetMultizoneIFaces();
	
	// Get the total number faces in the i--direction.
  const auto nFacesIDir = idir_faces_multizone.GetnFacesTotal(); 

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesIDir; i++)
	{
		// Get the relevant flattened index of the current face.
    //const auto face_info = multizone_geometry_container->GetFlattenedIndexIFace(i);
    const auto face_info = multizone_geometry_container->GetFlattenedIndexIFaceLoadBalanced(i);
    // Extract the current type of face.
    const auto face_type = face_info.mFaceType;
		
    // Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
        // Extract the current face.
        const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );

				const auto iZone  = family_face.GetiZone();
				const auto iElemL = internal_face.mIndexElementM;

        auto* physical_element = multizone_solver_container->GetSinglezoneSolver(iZone)->GetPhysicalElement(iElemL);

        auto& tmpL = physical_element->mResMinus;
        auto& resL = physical_element->mRes2D;
        
        for(size_t l=0; l<resL.size(); l++) resL[l] += tmpL[l];

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the family owning the current face.
        const auto& family_face = idir_faces_multizone.GetInterfaceFacesFamily( face_info.mIndexFamily );

        // Extract the actual face.
        const auto& interface_face = family_face.GetInterfaceFace( face_info.mIndexFace ); 

        const auto iZoneM = family_face.GetiZone();
        const auto iZoneP = family_face.GetjZone();

        const auto iElemM = interface_face.mIndexElementI;
        const auto iElemP = interface_face.mIndexElementJ;

        auto* physical_element_m = multizone_solver_container->GetSinglezoneSolver(iZoneM)->GetPhysicalElement(iElemM);
        auto* physical_element_p = multizone_solver_container->GetSinglezoneSolver(iZoneP)->GetPhysicalElement(iElemP);

        auto& tmp_m = physical_element_m->mResMinus;
        auto& res_m = physical_element_m->mRes2D;

        auto& tmp_p = physical_element_p->mResMinus;
        auto& res_p = physical_element_p->mRes2D;

        if( family_face.GetisMaxFaceI() )
        {
          for(size_t l=0; l<res_m.size(); l++) res_m[l] += tmp_m[l];
        }

        if( family_face.GetisMaxFaceJ() )
        {
          for(size_t l=0; l<res_p.size(); l++) res_p[l] += tmp_p[l];
        }

        break;
      }

			default: ERROR("Unknown IFace type, could not accumulate its residual.");
		}
	}
}

//-----------------------------------------------------------------------------------

void NResidualComputation::AccumulateJFaceResidualsCollective
(
 const CMultizoneGeometry *multizone_geometry_container,
 CMultizoneSolver         *multizone_solver_container
)
 /*
  *
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

	// Get the total faces in the j-direction.
	const auto& jdir_faces_multizone = multizone_geometry_container->GetMultizoneJFaces();
	
	// Get the total number faces in the  j-direction.
  const auto nFacesJDir = jdir_faces_multizone.GetnFacesTotal(); 

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nFacesJDir; i++)
	{
		// Get the relevant flattened index of the current face.
    //const auto face_info = multizone_geometry_container->GetFlattenedIndexJFace(i);
    const auto face_info = multizone_geometry_container->GetFlattenedIndexJFaceLoadBalanced(i);
    // Extract the current type of face.
    const auto face_type = face_info.mFaceType;

		// Check what type of face we are dealing with.
		switch( face_type )
		{
			case( ETypeFaceGeometry::INTERNAL ):
			{
        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInternalFacesFamily( face_info.mIndexFamily );
       
        // Extract the current face.
        const auto& internal_face = family_face.GetInternalFace( face_info.mIndexFace );

				const auto iZone  = family_face.GetiZone();
				const auto iElemB = internal_face.mIndexElementM;

        auto* physical_element = multizone_solver_container->GetSinglezoneSolver(iZone)->GetPhysicalElement(iElemB);

        auto& tmpB = physical_element->mResMinus;
        auto& resB = physical_element->mRes2D;
        
        for(size_t l=0; l<resB.size(); l++) resB[l] += tmpB[l];

				break;
			}

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the family owning the current face.
        const auto& family_face = jdir_faces_multizone.GetInterfaceFacesFamily( face_info.mIndexFamily );

        // Extract the actual face.
        const auto& interface_face = family_face.GetInterfaceFace( face_info.mIndexFace ); 

        const auto iZoneM = family_face.GetiZone();
        const auto iZoneP = family_face.GetjZone();

        const auto iElemM = interface_face.mIndexElementI;
        const auto iElemP = interface_face.mIndexElementJ;

        auto* physical_element_m = multizone_solver_container->GetSinglezoneSolver(iZoneM)->GetPhysicalElement(iElemM);
        auto* physical_element_p = multizone_solver_container->GetSinglezoneSolver(iZoneP)->GetPhysicalElement(iElemP);

        auto& tmp_m = physical_element_m->mResMinus;
        auto& res_m = physical_element_m->mRes2D;

        auto& tmp_p = physical_element_p->mResMinus;
        auto& res_p = physical_element_p->mRes2D;

        if( family_face.GetisMaxFaceI() )
        {
          for(size_t l=0; l<res_m.size(); l++) res_m[l] += tmp_m[l];
        }

        if( family_face.GetisMaxFaceJ() )
        {
          for(size_t l=0; l<res_p.size(); l++) res_p[l] += tmp_p[l];
        }

        break;
      }

			default: ERROR("Unknown JFace type, could not update its residual.");
		}
	}
}

//-----------------------------------------------------------------------------------

void NResidualComputation::ApplyInverseMassMatricesCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray
)
 /*
  *
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

	// Get the total number of elements in all zones.
	const auto nElemTotal = multizone_geometry_container->GetnElemTotal();

  // Borrow memory once, and use all the existing work array (which is more than enough).
  auto tmp = workarray.GetWorkMatrixAS3( 1, workarray.size() );

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
    // Extract the element indices.
		const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);

		// Deduce the current element's zone and index.
		const auto iZone = elem_info.mIndexZone; 
		const auto iElem = elem_info.mIndexElem; 

		// Extract the relevant solver.
		auto* solver  = multizone_solver_container->GetSinglezoneSolver(iZone);
		// Extract the relevant physical element.
		auto* element = solver->GetPhysicalElement(iElem);

    const as3double *minv = element->mInvMassMatrix.data();
    as3double        *res = element->mRes2D.data();

    solver->GetTensorProduct()->CompileTimeApplyInverseMassMatrix(minv, res, tmp.data()); 
  }
}



