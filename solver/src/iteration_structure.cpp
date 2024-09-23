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
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3vector1d<as3double>                   &monitordata,
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
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3vector1d<as3double>                   &monitordata,
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
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 CPoolMatrixAS3<as3double>                &workarray,
 as3vector1d<as3double>                   &monitordata,
 as3double                                 localtime
)
 /*
	* Function that computes the residual in all zones.
	*/
{
	// Extract the number of zones.
	const unsigned short nZone = config_container->GetnZone();

	// For now, loop over the zones and naively update the residuals.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		// Extract pointer to the relevant solver.
		auto* solver = solver_container[iZone].get();
		

		// Compute all volume terms. Note, this step also initializes the residual.
		solver->ComputeVolumeResidual(localtime, monitordata, workarray);

		// Compute all surface terms in the i-direction.
		solver->ComputeSurfaceResidualIDir(solver_container, geometry_container, 
				                               monitordata, workarray, localtime);

		// Compute all surface terms in the j-direction.
		solver->ComputeSurfaceResidualJDir(solver_container, geometry_container, 
				                               monitordata, workarray, localtime);

	}

	// TODO
	// Loop over the interface boundaries, if need be.
	for( auto& interface: interface_container )
	{
		interface->ComputeInterfaceResidual(solver_container, workarray);
	}
}

//-----------------------------------------------------------------------------------

void CIteration::GridSweep
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container, 
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 as3vector1d<as3double>                   &monitordata,
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
											solver_container,
											interface_container,
											workarray,
											monitordata,
											localtime);


	// Compute the residual over all zones.
	ComputeResidual(config_container,
			            geometry_container,
									solver_container,
									interface_container,
									workarray,
									monitordata,
									localtime);


	// Check for any postprocessing steps.
	PostProcessIteration(config_container,
			                 geometry_container,
											 solver_container,
											 interface_container,
											 workarray,
											 monitordata,
											 localtime);
}
