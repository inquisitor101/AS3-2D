#include "temporal_structure.hpp"

#include <cstddef>


//-----------------------------------------------------------------------------------
// ITemporal member functions.
//-----------------------------------------------------------------------------------


ITemporal::ITemporal
(
 const CConfig *config_container
)
 /*
	* Constructor for the temporal-method interface.
	*/
{

}

//-----------------------------------------------------------------------------------

ITemporal::~ITemporal
(
 void
)
 /*
	* Destructor for ITemporal.
	*/
{

}


//-----------------------------------------------------------------------------------
// CSSPRK3Temporal member functions.
//-----------------------------------------------------------------------------------


CSSPRK3Temporal::CSSPRK3Temporal
(
 const CConfig *config_container
)
	:
		ITemporal(config_container)
 /*
	* Constructor that sets the SSP-RK3 stage coefficients.
	*/
{
	// Initialize the SSP-RK3 coefficients. 
	mRk3a[0] = static_cast<as3double>(1.0);
	mRk3a[1] = static_cast<as3double>(3.0/4.0);
	mRk3a[2] = static_cast<as3double>(1.0/3.0);

	mRk3b[0] = static_cast<as3double>(1.0);
	mRk3b[1] = static_cast<as3double>(1.0/4.0);
	mRk3b[2] = static_cast<as3double>(2.0/3.0);

	mRk3c[0] = static_cast<as3double>(0.0);
	mRk3c[1] = static_cast<as3double>(1.0);
	mRk3c[2] = static_cast<as3double>(1.0/2.0);
}

//-----------------------------------------------------------------------------------

CSSPRK3Temporal::~CSSPRK3Temporal
(
 void
)
 /*
	* Destructor for CSSPRK3Temporal.
	*/
{

}

//-----------------------------------------------------------------------------------

void CSSPRK3Temporal::UpdateTime
(
 CConfig            *config_container,
 CMultizoneGeometry *multizone_geometry_container,
 CIteration         *iteration_container,
 COpenMP            *openmp_container,
 CMultizoneSolver   *multizone_solver_container,
 as3double           physicaltime,
 as3double           timestep
)
 /*
	* Function that advances one time step through the three SSP-RK3 stages.
	*/
{
	// Loop over all Runge-Kutta stages.
	for(std::size_t iStageRK=0; iStageRK<mNStageRK; iStageRK++)
	{
		// Local physical time.
		const as3double localtime = physicaltime + mRk3c[iStageRK]*timestep;

    // The iteration owns traversal, workspace, and synchronization for this stage.
    iteration_container->AdvanceStage(multizone_geometry_container,
                                      multizone_solver_container,
                                      this, localtime, timestep, iStageRK);
	}
}

//-----------------------------------------------------------------------------------

void CSSPRK3Temporal::UpdateElementSolution
(
 CPhysicalElement *physical_element,
 as3double         timestep,
 std::size_t       iStageRK
) const
 /*
	* Function that applies one SSP-RK3 stage formula to one element.
	*/
{
	// Local RK coefficients.
	const as3double alpha = mRk3a[iStageRK];
	const as3double beta  = mRk3b[iStageRK];
	const as3double bdt   = beta * timestep;

	// Extract the solution, spatial residual, and temporal history.
	auto& sol = physical_element->mSol2D;
	auto& res = physical_element->mRes2D;
	auto& old = physical_element->mSolOld2D;

	// Preserve the first-stage distinction and the original update expressions.
	if( iStageRK == 0 )
	{
#pragma omp simd
		for(std::size_t l=0; l<sol.size(); l++)
		{
			old[l]  = sol[l];
			sol[l] += bdt*res[l];
		}
	}
	else
	{
		const as3double oma = C_ONE - alpha;
#pragma omp simd
    for(std::size_t l=0; l<sol.size(); l++)
			sol[l] = oma*sol[l] + alpha*old[l] + bdt*res[l];
	}
}
