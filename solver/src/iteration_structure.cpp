#include "iteration_structure.hpp"


//-----------------------------------------------------------------------------------
// CIteration member functions.
//-----------------------------------------------------------------------------------


CIteration::CIteration
(
 const CConfig          *config_container,
 const CMultizoneSolver *multizone_solver_container
)
 /*
	* Constructor for the iteration class.
	*/
{
	// Check the relevant number of data entries required in the work array.
	for(unsigned short iZone=0; iZone<config_container->GetnZone(); iZone++)
	{
    // Get reference to the relevant solver.
    const auto* solver_container = multizone_solver_container->GetSinglezoneSolver(iZone);

		// For now, take the maximum number of items.
		switch( config_container->GetTypeSolver(iZone) )
		{
			case(ETypeSolver::EE):
			{
				const size_t nItem2D = 3;  // volume  terms needed.
				const size_t nItem1D = 2;  // surface terms needed.
				
				const size_t nVar    = solver_container->GetnVar();
				const size_t nInt1D  = solver_container->GetStandardElement()->GetnInt1D();
				const size_t nInt2D  = solver_container->GetStandardElement()->GetnInt2D();

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
 CConfig                   *config_container,
 CGeometry                 *geometry_container,
 COpenMP                   *openmp_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime 
)
 /*
	* Function that preprocesses the solution, before sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::PostProcessIteration
(
 CConfig                   *config_container,
 CGeometry                 *geometry_container,
 COpenMP                   *openmp_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime 
)
 /*
	* Function that postprocesses the solution, after sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::ComputeResiduals
(
 CConfig                   *config_container,
 CGeometry                 *geometry_container,
 COpenMP                   *openmp_container,
 CMultizoneSolver          *multizone_solver_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
)
 /*
	* Function that computes the residual in all zones. Note, these all use 
  * a collective call via an optimized OpenMP for loop. They must be executed 
  * inside an OpenMP parallel region
	*/
{
  // First, we compute the volume residual, which also initializes the residuals.
  NResidualComputation::ComputeVolumeResidualsCollective(geometry_container, 
                                                         multizone_solver_container, 
                                                         workarray, localtime);

  // Then, we compute the IFace residuals.
  NResidualComputation::ComputeIFaceResidualsCollective(geometry_container, 
                                                        multizone_solver_container, 
                                                        workarray, localtime);

  // Afterwards, we must accumulate the temporary stored IFace residuals.
  NResidualComputation::AccumulateIFaceResidualsCollective(geometry_container, 
                                                           multizone_solver_container);

  // Same with JFace residuals.
  NResidualComputation::ComputeJFaceResidualsCollective(geometry_container,
                                                        multizone_solver_container,
                                                        workarray, localtime);

  // Also, accumulate the temporary stored JFace residuals.
  NResidualComputation::AccumulateJFaceResidualsCollective(geometry_container,
                                                           multizone_solver_container);

  // Finally, we update the residuals by including the mass matrix's effect.
  NResidualComputation::ApplyInverseMassMatricesCollective(geometry_container, 
                                                           multizone_solver_container, 
                                                           workarray);
}

//-----------------------------------------------------------------------------------

void CIteration::GridSweep
(
 CConfig          *config_container,
 CGeometry        *geometry_container,
 COpenMP          *openmp_container,
 CMultizoneSolver *multizone_solver_container, 
 as3double         localtime 
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
											multizone_solver_container,
											workarray,
											localtime);


	// Compute the residual over all zones.
	ComputeResiduals(config_container,
			             geometry_container,
									 openmp_container,
									 multizone_solver_container,
									 workarray,
									 localtime);


	// Check for any postprocessing steps.
	PostProcessIteration(config_container,
			                 geometry_container,
											 openmp_container,
											 multizone_solver_container,
											 workarray,
											 localtime);
}
