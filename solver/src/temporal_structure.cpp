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


