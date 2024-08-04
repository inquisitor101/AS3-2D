#include "marker_structure.hpp"


//-----------------------------------------------------------------------------------
// CMarker member functions.
//-----------------------------------------------------------------------------------


CMarker::CMarker
(
 unsigned short iZone,
 std::string    tagname,
 ETypeBCs       markertype
)
	:
		mZoneID(iZone),
		mNameMarkerTag(tagname),
		mTypeMarker(markertype)
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


//-----------------------------------------------------------------------------------
// CInterfaceMarker member functions.
//-----------------------------------------------------------------------------------


CInterfaceMarker::CInterfaceMarker
(
 unsigned short iZone,
 std::string    tagnameI,
 unsigned short jZone,
 std::string    tagnameJ,
 ETypeBCs       markertype
)
	:
		CMarker(iZone, tagnameI, markertype),
		mMatchingZoneID(jZone), mMatchingNameMarkerTag(tagnameJ)
 /*
	* Constructor for the interface-type marker.
	*/
{
	
}

//-----------------------------------------------------------------------------------

CInterfaceMarker::~CInterfaceMarker
(
 void
)
 /*
	* Destructor, which cleans up after the driver class.
	*/
{

}


