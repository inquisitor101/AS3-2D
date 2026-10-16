#include "boundary_structure.hpp"


//-----------------------------------------------------------------------------------
// IBoundary member functions.
//-----------------------------------------------------------------------------------


IBoundary::IBoundary
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 unsigned int              index
)
	:
		mIndex(index)
 /*
	* Constructor for the interface boundary class.
	*/
{

}

//-----------------------------------------------------------------------------------

IBoundary::~IBoundary
(
 void
)
 /*
	* Destructor, which cleans up after the interface boundary class.
	*/
{

}





