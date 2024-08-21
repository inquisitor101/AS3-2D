#include "marker_structure.hpp"


//-----------------------------------------------------------------------------------
// CMarker member functions.
//-----------------------------------------------------------------------------------


CMarker::CMarker
(
 unsigned short            iZone,
 std::string               tagname,
 ETypeBC                   typemarker,
 as3vector1d<unsigned int> elements,
 as3vector1d<EFaceElement> faces 
)
	:
		mZoneID(iZone),
		mNameMarkerTag(tagname),
	  mMarkerBC(typemarker),
		mElementIndices(elements),
		mElementFaces(faces)	
 /*
	* Constructor for the marker geometry, which contains a single marker region.
	*/
{
	
}

//-----------------------------------------------------------------------------------

CMarker::~CMarker
(
 void
)
 /*
	* Destructor, which cleans up after the driver class.
	*/
{

}


