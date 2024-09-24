#include "temporal_structure.hpp"


//-----------------------------------------------------------------------------------
// ITemporal member functions.
//-----------------------------------------------------------------------------------


ITemporal::ITemporal
(
 CConfig *config_container
)
 /*
	* Constructor for the interface temporal class.
	*/
{

}

//-----------------------------------------------------------------------------------

ITemporal::~ITemporal
(
 void
)
 /*
	* Destructor, which cleans up after the interface temporal class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CSSPRK3Temporal member functions.
//-----------------------------------------------------------------------------------


CSSPRK3Temporal::CSSPRK3Temporal
(
 CConfig *config_container
)
	:
		ITemporal(config_container)
 /*
	* Constructor for the strong-stability-preserving 3rd-order Runge-Kutta (SSP-RK3) temporal class.
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
	* Destructor, which cleans up after the SSP-RK3 temporal class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CSSPRK3Temporal::UpdateTime
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 CIteration                               *iteration_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 as3double                                 physicaltime,
 as3double                                 timestep,
 as3vector1d<as3double>                   &monitordata
)
 /*
	* Function that computes the upcoming solution in time, based on a SSP-RK3.
	*/
{
	// Loop over all Runge-Kutta stages.
	for(unsigned short iStageRK=0; iStageRK<mNStageRK; iStageRK++)
	{
		// Local physical time.
		const as3double localtime = physicaltime + mRk3c[iStageRK]*timestep;

    // Perform a single stage evaluation, based on a SSP-RK3 .
    EvaluateSSPRK3(config_container,
				           geometry_container,
               		 iteration_container,
               		 solver_container,
									 interface_container,
               		 localtime, timestep,
									 iStageRK,
               		 monitordata);
	}
}

//-----------------------------------------------------------------------------------

void CSSPRK3Temporal::EvaluateSSPRK3
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 CIteration                               *iteration_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container,
 as3double                                 localtime,
 as3double                                 timestep,
 unsigned short                            iStageRK,
 as3vector1d<as3double>                   &monitordata
)
 /*
	* Function that does a single stage SSP-RK3 evaluation.
	*/
{
	// Local RK coefficients.
	const as3double alpha = mRk3a[iStageRK];
	const as3double beta  = mRk3b[iStageRK];

	// First, compute the residual.
	iteration_container->GridSweep(config_container,
			                           geometry_container,
																 solver_container,
																 interface_container,
																 monitordata,
																 localtime);


	// Update the solution in time. For computational efficiency, make a distinction 
	// between the first stage (iStage=0) and the remaining ones.
	if( iStageRK == 0 )
	{
		// Abbreviation.
		const as3double bdt = beta*timestep;

		// Loop over each solver in every zone.
		for( auto& solver: solver_container )
		{
			// Loop over every element in this zone.
			for( auto& element: solver->GetPhysicalElement() )
			{
				// Extract current solution and residual.
				auto& sol = element->mSol2D;
				auto& res = element->mRes2D;
		
				// Extract the tentative (sol) solution state.
				auto& old = element->mSolOld2D;

				// Loop over each node and compute its updated value.
				for(size_t l=0; l<sol.size(); l++)
				{
					old[l]  = sol[l];
					sol[l] += bdt*res[l];
				}
			}
		}
	}
	else
	{
		// Abbreviation.
		const as3double bdt = beta*timestep;
		const as3double oma = C_ONE - alpha;

		// Loop over each solver in every zone.
		for( auto& solver: solver_container )
		{
			// Loop over every element in this zone.
			for( auto& element: solver->GetPhysicalElement() )
			{
				// Extract current solution and residual.
				auto& sol = element->mSol2D;
				auto& res = element->mRes2D;
		
				// Extract the tentative (sol) solution state.
				auto& old = element->mSolOld2D;

				// Loop over each node and compute its updated value.
				for(size_t l=0; l<sol.size(); l++)
				{
					sol[l] = oma*sol[l] + alpha*old[l] + bdt*res[l];
				}
			}
		}
	}
}




