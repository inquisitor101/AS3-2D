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
	// Number of stages is three.
	const size_t nStageRK = 3;

	// Initialize the SSP-RK3 coefficients.
	rk3a.resize(nStageRK);
	rk3b.resize(nStageRK);
	rk3c.resize(nStageRK);

	// Initialize the SSP-RK3 coefficients. 
	rk3a[0] = 1.0;
	rk3a[1] = 3.0/4.0;
	rk3a[2] = 1.0/3.0;

	rk3b[0] = 1.0;
	rk3b[1] = 1.0/4.0;
	rk3b[2] = 2.0/3.0;

	rk3c[0] = 0.0;
	rk3c[1] = 1.0;
	rk3c[2] = 1.0/2.0;
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
	// Extract number of RK stages in a SSP-RK3.
	const size_t nStageRK = rk3a.size();

	// Loop over all Runge-Kutta stages.
	for(size_t iStageRK=0; iStageRK<nStageRK; iStageRK++)
	{
		// Local physical time.
		const as3double localtime = physicaltime + rk3c[iStageRK]*timestep;

		// Local RK coefficients.
		const as3double alpha = rk3a[iStageRK];
		const as3double beta  = rk3b[iStageRK];

    // Perform a single stage evaluation, based on a SSP-RK3 .
    EvaluateSSPRK3(config_container,
				           geometry_container,
               		 iteration_container,
               		 solver_container,
									 interface_container,
               		 localtime, timestep,
               		 alpha, beta,
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
 as3double                                 alpha,
 as3double                                 beta,
 as3vector1d<as3double>                   &monitordata
)
 /*
	* Function that does a single stage SSP-RK3 evaluation.
	*/
{
	// First, compute the residual.
	iteration_container->GridSweep(config_container,
			                           geometry_container,
																 solver_container,
																 interface_container,
																 monitordata,
																 localtime);
}




