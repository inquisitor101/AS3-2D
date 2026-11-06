#include "geometry/CElementGeometry.hpp"

#include "common/AS3Constants.hpp"

#include <algorithm>
#include <cmath>


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

//-----------------------------------------------------------------------------------

bool CElementGeometry::CoordinatesMatch
(
 as3double xi,
 as3double xj,
 as3double translation
)
 /*
	* Function that allows for roundoff in the coordinates and their translation.
	*/
{
	constexpr as3double atol           = static_cast<as3double>(1.0e-10);
	constexpr as3double rtol           = static_cast<as3double>(1.0e-8);
	constexpr as3double roundoffFactor = static_cast<as3double>(4)*C_EPS;

	const as3double translated = xi + translation;
	const as3double difference = std::abs(translated - xj);
	const as3double tolerance  = atol + rtol*std::max(std::abs(translated), std::abs(xj));

	// Use the original coordinate and displacement scales in case the sum cancels.
	const as3double scale    = std::max({std::abs(xi), std::abs(xj), std::abs(translation)});
	const as3double roundoff = roundoffFactor*scale;
	return std::isfinite(difference) && difference <= std::max(tolerance, roundoff);
}
