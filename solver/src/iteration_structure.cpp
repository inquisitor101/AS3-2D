#include "iteration_structure.hpp"

#include <algorithm>
#include <cstddef>
#include "temporal_structure.hpp"


//-----------------------------------------------------------------------------------
// CIteration member functions.
//-----------------------------------------------------------------------------------


CIteration::CIteration
(
 const CConfig          *config_container,
 const CMultizoneSolver *multizone_solver_container
)
 /*
	* Constructor that sizes the workspace needed by each solver thread.
	*/
{
	// Check the relevant number of data entries required in the work array.
	for(std::size_t iZone=0; iZone<config_container->GetnZone(); iZone++)
	{
    // Get reference to the relevant solver.
    const auto* solver_container = multizone_solver_container->GetSinglezoneSolver(iZone);

		// For now, take the maximum number of items.
		switch( config_container->GetTypeSolver(iZone) )
		{
			case(ETypeSolver::EE):
			{
				const std::size_t nItem2D = 3;  // volume  terms needed.
				const std::size_t nItem1D = 2;  // surface terms needed.
				
				const std::size_t nVar    = solver_container->GetnVar();
				const std::size_t nInt1D  = solver_container->GetStandardElement()->GetnInt1D();
				const std::size_t nInt2D  = solver_container->GetStandardElement()->GetnInt2D();

				// Compute the total number of required volume and surface terms in the work array.
				const std::size_t nVol  = nItem2D*nInt2D*nVar;
				const std::size_t nSurf = nItem1D*nInt1D*nVar;

				// Take the maximum storage between the volume and surface terms.
				const std::size_t nData = std::max( nVol, nSurf );

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
	* Destructor for CIteration.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::ComputeResiduals
(
 CConfig                   *config_container,
 CMultizoneGeometry        *multizone_geometry_container,
 COpenMP                   *openmp_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
	* Function that collectively computes surface and element residuals for all zones.
	*/
{
  // First, compute the surface residuals into their local face buffers.
  NResidualComputation::ComputeSurfaceResidualsCollectiveNowait(multizone_geometry_container,
                                                               multizone_solver_container,
                                                               workarray, localtime);

  // Then, compute the volume terms, gather the surfaces, and apply the inverse mass matrix per element.
  NResidualComputation::ComputeElementResidualsCollective(multizone_geometry_container,
                                                          multizone_solver_container,
                                                          workarray, localtime);
}

//-----------------------------------------------------------------------------------

void CIteration::AdvanceStage
(
 const CMultizoneGeometry *multizone_geometry_container,
 CMultizoneSolver         *multizone_solver_container,
 const ITemporal          *temporal_container,
 as3double                 localtime,
 as3double                 timestep,
 std::size_t            iStageRK
)
 /*
  * Function that computes surface residuals, synchronizes, then assembles and updates each element.
  */
{
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Each thread owns its workspace for both surface and volume computations.
    CPoolMatrixAS3<as3double> workarray(mNWorkData);

    // Compute the surface residuals into their local face buffers.
    NResidualComputation::ComputeSurfaceResidualsCollectiveNowait(multizone_geometry_container,
                                                                 multizone_solver_container,
                                                                 workarray, localtime);

    // Wait for all surface contributions before changing any element solution.
#ifdef HAVE_OPENMP
#pragma omp barrier
#endif

    // Get the total number of elements in all zones.
    const auto nElemTotal = multizone_geometry_container->GetnElemTotal();

#ifdef HAVE_OPENMP
#pragma omp for schedule(static) nowait
#endif
    for(std::size_t i=0; i<nElemTotal; i++)
    {
      // Extract the element indices.
      const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);
      const auto iZone = elem_info.mIndexZone;
      const auto iElem = elem_info.mIndexElem;

      // Extract the relevant solver, grid, and physical element.
      auto* solver  = multizone_solver_container->GetSinglezoneSolver(iZone);
      auto* grid    = multizone_geometry_container->GetSinglezoneGeometry(iZone);
      auto* element = solver->GetPhysicalElement(iElem);

      // Compute the volume terms, accumulate the surfaces, and apply the inverse mass matrix.
      NResidualComputation::ComputeElementResidual(grid, solver, workarray, localtime, iElem);

      // The temporal method supplies the update formula for this element.
      temporal_container->UpdateElementSolution(element, timestep, iStageRK);
    }
#ifdef ENABLE_NAN_CHECK
    // Floating-point status is thread-local; checking only the driver thread
    // would miss exceptions raised by the other OpenMP workers.
    NError::CheckFloatingError();
#endif
  } // The parallel-region barrier completes all solution updates before the next stage.
}

//-----------------------------------------------------------------------------------

void CIteration::GridSweep
(
 CConfig            *config_container,
 CMultizoneGeometry *multizone_geometry_container,
 COpenMP            *openmp_container,
 CMultizoneSolver   *multizone_solver_container,
 as3double           localtime
)
 /*
	* Function that collectively computes spatial residuals without a temporal update.
	*/
{
	// Initialize a separate work array for each thread in the enclosing parallel region.
	CPoolMatrixAS3<as3double> workarray(mNWorkData);

	// Compute the residual over all zones.
	ComputeResiduals(config_container,
			             multizone_geometry_container,
									 openmp_container,
									 multizone_solver_container,
									 workarray,
									 localtime);
#ifdef ENABLE_NAN_CHECK
	// Check floating-point errors raised by the current worker during the grid sweep.
	NError::CheckFloatingError();
#endif
}
