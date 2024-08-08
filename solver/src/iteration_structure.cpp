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
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
 as3double                              localtime, 
 CPoolMatrixAS3<as3double>             &workarray,
 as3vector1d<as3double>                &monitordata

)
 /*
	* Function that preprocesses the solution, before sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::PostProcessIteration
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
 as3double                              localtime, 
 CPoolMatrixAS3<as3double>             &workarray,
 as3vector1d<as3double>                &monitordata
)
 /*
	* Function that postprocesses the solution, after sweeping the grid.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIteration::ComputeResidual
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
 as3double                              localtime,
 CPoolMatrixAS3<as3double>             &workarray,
 as3vector1d<as3double>                &monitordata
)
 /*
	* Function that computes the residual in all zones.
	*/
{
	// NOTE
	// This is the most performance-dependent function, at least under the current design choices.
	// What we know:
	// 1) The volume terms are the most expensive to compute. However, they are compact terms, completely decoupled.
	// 2) The surface terms are coupled, as they rely on neighboring elements and boundary values.
	//
	// Probably, better to seperate the surface terms and the volume terms. 
	// Also, perhaps it is more efficient to consider internal and boundary surface faces separately?



	// Extract the number of zones.
	const unsigned short nZone = config_container->GetnZone();

	// For now, loop over the zones and naively update the residuals.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		// Extract pointer to the relevant solver.
		auto* solver       = solver_container[iZone].get();
		// Extract pointer to the relevant grid zone.
		auto* grid         = geometry_container->GetZoneGeometry(iZone);
		// Extract the number of elements in this zone.
		const size_t nElem = grid->GetnElem();

	

		// Loop over the elements and compute the residuals.
		for(size_t iElem=0; iElem<nElem; iElem++)
		{
			// Compute volume terms.
			solver->ComputeVolumeResidual(iElem, localtime, monitordata, workarray);

			// Compute surface terms.

			// Multiply by the inverse mass matrix.
		}
	}
}

//-----------------------------------------------------------------------------------

void CIteration::GridSweep
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
 as3double                              localtime, 
 as3vector1d<as3double>                &monitordata
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
											localtime,
											workarray,
											monitordata);


	// Compute the residual over all zones.
	ComputeResidual(config_container,
			            geometry_container,
									solver_container,
									localtime,
									workarray,
									monitordata);


	// Check for any postprocessing steps.
	PostProcessIteration(config_container,
			                 geometry_container,
											 solver_container,
											 localtime,
											 workarray,
											 monitordata);
}
