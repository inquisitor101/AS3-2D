#include "boundary_structure.hpp"


//-----------------------------------------------------------------------------------
// IBoundary member functions.
//-----------------------------------------------------------------------------------


IBoundary::IBoundary
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CMarker   *marker_container
)
	:
		mMarker( marker_container )
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


//-----------------------------------------------------------------------------------
// CPeriodicBoundary member functions.
//-----------------------------------------------------------------------------------


CPeriodicBoundary::CPeriodicBoundary
(
 CConfig   *config_container,
 CGeometry *geometry_container,
 CMarker   *marker_container
)
	:
		IBoundary(config_container, geometry_container, marker_container)
 /*
	* Constructor for the periodic boundary condition class.
	*/
{
	// Determine where the matching marker is located.
	FindMatchingMarker(config_container, geometry_container);

	// Ensure both markers are specified as periodic BCs.
	if( mMarker->GetMarkerBC()         != ETypeBC::PERIODIC ) ERROR("Parent marker must use periodic BCs.");
	if( mMatchingMarker->GetMarkerBC() != ETypeBC::PERIODIC ) ERROR("Matching marker must use periodic BCs."); 
}

//-----------------------------------------------------------------------------------

CPeriodicBoundary::~CPeriodicBoundary
(
 void
)
 /*
	* Destructor, which cleans up after the periodic boundary condition class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CPeriodicBoundary::FindMatchingMarker
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that finds the matching marker for this periodic boundary.
	*/
{
	// Extract the name of the marker we are looking for.
	auto jmark = config_container->FindMatchingPeriodicMarker( mMarker->GetNameMarkerTag() ); 

	// Flag on whether the marker has been found or not.
	bool found = false;

	// Loop over the possible markers and find the matching marker.
	for(unsigned short iZone=0; iZone<config_container->GetnZone(); iZone++)
	{
		// Extract current zone geometry.
		auto* zone = geometry_container->GetZoneGeometry(iZone);

		// Loop over each marker and check for the marker name.
		for( auto& m: zone->GetMarker() )
		{
			if( m->GetNameMarkerTag() == jmark ) 
			{
				// Set the pointer to the matching marker and break from the loop.
				mMatchingMarker = m.get(); found = true; break;
			}
		}
		if( found ) break;
	}

	// Issue an error, if the matching marker couldn't be located.
	if( !found ) ERROR("Could not be found the matching periodic marker: " + jmark);
}


