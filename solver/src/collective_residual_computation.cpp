#include "collective_residual_computation.hpp"

#include <cstddef>


//-----------------------------------------------------------------------------------
// Local surface residual functions.
//-----------------------------------------------------------------------------------

static void ComputeIFaceResidualsCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  * Function that collectively computes I-face residuals without an end barrier.
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  NResidualComputation::AssertParallelRegion();

	// Get the total faces in the i-direction.
	const auto& idir_faces_multizone = multizone_geometry_container->GetMultizoneIFaces();
	
	// Get the total number faces in the i--direction.
  const auto nFacesIDir = idir_faces_multizone.GetnFacesTotal(); 

#ifdef HAVE_OPENMP
#pragma omp for schedule(static) nowait
#endif
	for(std::size_t i=0; i<nFacesIDir; i++)
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
        // Extract the batch owning the current face.
        const auto& face_batch = idir_faces_multizone.GetInteriorFaces().GetBatch( face_info.mIndexBatch );
        // Extract the current face.
				const auto& internal_face = face_batch.GetFace( face_info.mIndexFace );
			
        const auto iZone  = face_batch.GetZoneIndex();
				const auto iElemL = internal_face.mIndexElementMinus;

				const auto* grid = multizone_geometry_container->GetSinglezoneGeometry(iZone);

				// Compute the internal face residual.
				multizone_solver_container->GetSinglezoneSolver(iZone)->ComputeSurfaceResidualIDir(grid, workarray, localtime, iElemL);

				break;
			}

      case( ETypeFaceGeometry::BOUNDARY ):
      {
        // Extract the relevant batch index.
        const auto iBatch = face_info.mIndexBatch;
        // Extract the associated boundary index.
        const auto iBoundary = multizone_solver_container->GetIndexBoundaryIFace(iBatch);

        // Extract the associated boundary object.
        const auto* boundary = multizone_solver_container->GetBoundary(iBoundary);
        // Extract the batch owning the current face.
        const auto& face_batch = idir_faces_multizone.GetBoundaryFaces().GetBatch( iBatch );

        // Extract the relevant zone index.
        const auto iZone = face_batch.GetZoneIndex();
        
        // Extract the relevant grid and solver.
        const auto* grid = multizone_geometry_container->GetSinglezoneGeometry(iZone);
        auto* solver     = multizone_solver_container->GetSinglezoneSolver(iZone);

        // Compute the boundary face residual.
        boundary->ComputeBoundaryResidual(grid, solver, face_batch, face_info, workarray, localtime);       

        break;
      }

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the relevant batch index, which is the same as the interface container's index.
        const auto iBatch = face_info.mIndexBatch;
        // Extract the associated interface index.
        const auto iInterface = multizone_solver_container->GetIndexInterfaceIFace(iBatch); 
        
        // Extract the associated interface object.
        const auto* interface = multizone_solver_container->GetInterface(iInterface); 
        // Extract the batch owning the current face.
        const auto& face_batch = idir_faces_multizone.GetInterfaceFaces().GetBatch( iBatch );
       
        // Compute the interface residual.
        interface->ComputeInterfaceResidual(multizone_solver_container, face_batch, face_info, workarray, localtime);
        
        break;
      }

			default: ERROR("Unknown IFace type, could not compute its residual.");
		}
	}
}

//-----------------------------------------------------------------------------------

static void ComputeJFaceResidualsCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  * Function that collectively computes J-face residuals without an end barrier.
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  NResidualComputation::AssertParallelRegion();

	// Get the total faces in the j-direction.
	const auto& jdir_faces_multizone = multizone_geometry_container->GetMultizoneJFaces();
	
	// Get the total number faces in the  j-direction.
  const auto nFacesJDir = jdir_faces_multizone.GetnFacesTotal(); 

#ifdef HAVE_OPENMP
#pragma omp for schedule(static) nowait
#endif
	for(std::size_t i=0; i<nFacesJDir; i++)
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
        // Extract the batch owning the current face.
        const auto& face_batch = jdir_faces_multizone.GetInteriorFaces().GetBatch( face_info.mIndexBatch );
       
        // Extract the current face.
        const auto& internal_face = face_batch.GetFace( face_info.mIndexFace );

        const auto iZone  = face_batch.GetZoneIndex();
        const auto iElemB = internal_face.mIndexElementMinus;

        const auto* grid = multizone_geometry_container->GetSinglezoneGeometry(iZone);

				// Compute the internal face residual. 
				multizone_solver_container->GetSinglezoneSolver(iZone)->ComputeSurfaceResidualJDir(grid, workarray, localtime, iElemB);

				break;
			}

      case( ETypeFaceGeometry::BOUNDARY ):
      {
        // Extract the relevant batch index.
        const auto iBatch = face_info.mIndexBatch;
        // Extract the associated boundary index.
        const auto iBoundary = multizone_solver_container->GetIndexBoundaryJFace(iBatch);

        // Extract the associated boundary object.
        const auto* boundary = multizone_solver_container->GetBoundary(iBoundary);
        // Extract the batch owning the current face.
        const auto& face_batch = jdir_faces_multizone.GetBoundaryFaces().GetBatch( iBatch );
        // Extract the relevant zone index.
        const auto iZone = face_batch.GetZoneIndex();
        
        // Extract the relevant grid and solver.
        const auto* grid = multizone_geometry_container->GetSinglezoneGeometry(iZone);
        auto* solver     = multizone_solver_container->GetSinglezoneSolver(iZone);

        // Compute the boundary face residual.
        boundary->ComputeBoundaryResidual(grid, solver, face_batch, face_info, workarray, localtime);       

        break;
      }

      case( ETypeFaceGeometry::INTERFACE ):
      {
        // Extract the relevant batch index, which is the same as the interface container's index.
        const auto iBatch = face_info.mIndexBatch;
        // Extract the associated interface index.
        const auto iInterface = multizone_solver_container->GetIndexInterfaceJFace(iBatch); 
        
        // Extract the associated interface object.
        const auto* interface = multizone_solver_container->GetInterface(iInterface); 
        // Extract the batch owning the current face.
        const auto& face_batch = jdir_faces_multizone.GetInterfaceFaces().GetBatch( iBatch );

        // Compute the interface residual.
        interface->ComputeInterfaceResidual(multizone_solver_container, face_batch, face_info, workarray, localtime);
        
        break;
      }

			default: ERROR("Unknown JFace type, could not compute its residual.");
		}
	}
}

//-----------------------------------------------------------------------------------
// NResidualComputation namespace functions.
//-----------------------------------------------------------------------------------

void NResidualComputation::ComputeSurfaceResidualsCollectiveNowait
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  * Function that collectively computes both directional surface sweeps without an end barrier.
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

  // Compute both directions into their respective local face buffers.
  ComputeIFaceResidualsCollective(multizone_geometry_container,
                                  multizone_solver_container,
                                  workarray, localtime);

  ComputeJFaceResidualsCollective(multizone_geometry_container,
                                  multizone_solver_container,
                                  workarray, localtime);
}

//-----------------------------------------------------------------------------------

void NResidualComputation::ComputeElementResidualsCollective
(
 const CMultizoneGeometry  *multizone_geometry_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
  * Function that waits for surface writes and collectively assembles mass-scaled element residuals.
  */
{
  // Ensure this is an OpenMP parallel region, in case we're running in parallel.
  AssertParallelRegion();

  // Wait for all surface contributions before reading the local face buffers.
#ifdef HAVE_OPENMP
#pragma omp barrier
#endif

	// Get the total number of elements in all zones.
	const auto nElemTotal = multizone_geometry_container->GetnElemTotal();

#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
	for(std::size_t i=0; i<nElemTotal; i++)
	{
    // Extract the element indices.
		const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);

		// Deduce the current element's zone and index.
		const auto iZone = elem_info.mIndexZone; 
		const auto iElem = elem_info.mIndexElem; 

		// Extract the relevant solver.
		auto* solver  = multizone_solver_container->GetSinglezoneSolver(iZone);
		// Extract the relevant grid.
		auto* grid    = multizone_geometry_container->GetSinglezoneGeometry(iZone);

    // Compute the complete spatial residual for this element.
    ComputeElementResidual(grid, solver, workarray, localtime, iElem);
  }
}

//-----------------------------------------------------------------------------------

void NResidualComputation::ComputeElementResidual
(
 CSinglezoneGeometry       *singlezone_geometry_container,
 ISolver                   *singlezone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElem
)
 /*
  * Function that assembles one element residual and applies its inverse mass matrix.
  */
{
  // Extract the relevant physical element.
  auto* element = singlezone_solver_container->GetPhysicalElement(iElem);

  // Compute the volume terms first, which also initializes the element residual.
  singlezone_solver_container->ComputeVolumeResidual(singlezone_geometry_container, workarray, localtime, iElem);

  // Accumulate the temporary surface residuals into the element residual.
  auto& residual = element->mRes2D;
  for(const auto& surface : element->mResSurface)
    for(std::size_t l=0; l<residual.size(); l++) residual[l] += surface[l];

  // Borrow the workspace after the volume temporaries have been released.
  auto tmp = workarray.GetWorkMatrixAS3( 1, workarray.size() );

  const as3double *minv = element->mInvMassMatrix.data();
  as3double        *res = element->mRes2D.data();

  singlezone_solver_container->GetTensorProduct()->CompileTimeApplyInverseMassMatrix(minv, res, tmp.data());
}



