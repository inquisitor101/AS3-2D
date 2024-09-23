#include "marker_structure.hpp"


//-----------------------------------------------------------------------------------
// CMarker member functions.
//-----------------------------------------------------------------------------------


CMarker::CMarker
(
 unsigned short            zone,
 ETypeBC                   type,
 std::string               name,
 as3vector1d<EFaceElement> face,
 as3vector1d<unsigned int> mark
)
	:
		mZoneID(zone),
		mTypeBC(type),
		mNameMarker(name)
 /*
	* Constructor for the marker geometry, which contains a single marker region.
	*/
{
	// Ensure the faces and index markes are of the same size.
	if( face.size() != mark.size() ) ERROR("Mismatch in the size of face and element indices.");

	// Reserve memory for the element face properties.
	mElementFaces.reserve( face.size() );
	
	// Initialize the member variables.
	for(size_t i=0; i<face.size(); i++)
	{
		mElementFaces.emplace_back( mark[i], face[i] );
	}
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


