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
	* Constructor for the SSP-RK3 temporal class.
	*/
{

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


