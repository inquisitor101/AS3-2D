#include "geometry/CElementGeometry.hpp"


//-----------------------------------------------------------------------------------
// CElementGeometry member functions.
//-----------------------------------------------------------------------------------

CElementGeometry::CElementGeometry
(
 std::size_t nNode2D
)
  : mCoordSolDOFs(2, nNode2D)
 /*
  * Constructor that allocates geometry-node coordinates for one element.
  */
{

}

//-----------------------------------------------------------------------------------

CElementGeometry::~CElementGeometry
(
 void
)
 /*
	* Destructor for CElementGeometry.
	*/
{

}
