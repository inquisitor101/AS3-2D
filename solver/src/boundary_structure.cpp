#include "boundary_structure.hpp"


//-----------------------------------------------------------------------------------
// IBoundary member functions.
//-----------------------------------------------------------------------------------


IBoundary::IBoundary
(
 CConfig      *config_container,
 CGeometry    *geometry_container,
 CMarker      *marker_container,
 unsigned int  index
)
	:
		mIndex(index)
 /*
	* Constructor for the interface boundary class.
	*/
{
	// Ensure that the input marker is not null.
	if( !marker_container ) ERROR("Input marker_container cannot be nullptr.");
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





