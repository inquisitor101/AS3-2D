#include "spatial_structure.hpp"


//-----------------------------------------------------------------------------------
// ISpatial member functions.
//-----------------------------------------------------------------------------------


ISpatial::ISpatial
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor for the interface spatial discretization class.
	*/
{

}

//-----------------------------------------------------------------------------------

ISpatial::~ISpatial
(
 void
)
 /*
	* Destructor, which cleans up after the interface spatial discretization class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CEESpatial member functions.
//-----------------------------------------------------------------------------------


CEESpatial::CEESpatial
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
	:
		ISpatial(config_container, geometry_container)
 /*
	* Constructor for a pure Euler equations spatial class.
	*/
{

}

//-----------------------------------------------------------------------------------

CEESpatial::~CEESpatial
(
 void
)
 /*
	* Destructor, which cleans up after the Euler equations spatial class.
	*/
{

}


