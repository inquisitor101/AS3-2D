#include "face_structure.hpp"
#include "geometry_structure.hpp"


//-----------------------------------------------------------------------------------
// CZoneGeometry member functions.
//-----------------------------------------------------------------------------------


void CMultizoneFaceGeometry::InitializeInternalFaces
(
 const CGeometry *geometry_container
)
 /*
	*
	*/
{
	// Determine total number of internal faces in the i-direction.
	size_t nInternalFacesIDir = 0;
	for( const auto& zone : geometry_container->GetZoneGeometry() )
	{
		const size_t niElem = zone->GetnxElem();
    const size_t njElem = zone->GetnyElem();
	
		nInternalFacesIDir += njElem * (niElem-1);
	}
	
	// Reserve the needed amount of memory for the internal faces.
	mInternalFaces.reserve(nInternalFacesIDir);

	// Loop over the internal faces in all the zones and instantiate them.
  for( const auto& zone : geometry_container->GetZoneGeometry() )
  {
    const size_t  iZone = zone->GetZoneID();
    const size_t niElem = zone->GetnxElem();
    const size_t njElem = zone->GetnyElem();

    for(size_t j=0; j<njElem; j++)
    {
      for(size_t i=1; i<niElem; i++)
      {
        // Deduce the flattened element indices.
        const size_t ijElemR = j * niElem + i;
        const size_t ijElemL = ijElemR - 1;
        
				mInternalFaces.emplace_back( ijElemL, ijElemR, iZone ); 
      }
    }
	}

}

//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::InitializeInterfaceFaces
(
 const CConfig   *config_container,
 const CGeometry *geometry_container
)
 /*
	*
	*/
{
	// Extract the interface boundaries, specified by the user.
	auto& param_interface = config_container->GetInterfaceParamMarker();

	// If there are no interfaces specified, return without any initialization.
	if( param_interface.empty() ) return;

	// Allocate memory for the family of interfaces.
	mTotalInterfaceBoundary.reserve( param_interface.size() );

	// Temporary lambda to search for the zone of a given marker name.
	auto lFindMarker = [=](std::string &name) -> CMarker*
	{
		for( auto& zone: geometry_container->GetZoneGeometry() )
		{
			for( auto& marker: zone->GetMarker() )
			{
				if( marker->GetNameMarker() == name )
				{
					return marker.get();
				}
			}
		}
		
		// The program should not reach here, otherwise we have an error.
		ERROR("Could not find the interface marker.");

		// To avoid compiler problems, return something.
		return nullptr;
	};


	// Loop over all the zones and search for the matching markers and pair them.
	for( const auto& zone : geometry_container->GetZoneGeometry() )
	{

		// Get a pointer to the two markers forming this interface.
		const CMarker *imarker_container = lFindMarker(param_interface->mName);
		const CMarker *jmarker_container = lFindMarker(param_interface->mNameMatching);

		// Extract the zone ID of these markers.
		const unsigned short iZone = imarker_container->GetZoneID();
		const unsigned short jZone = jmarker_container->GetZoneID();


		// Temporary lambda to find the face location on a given marker.
		auto lFaceLocation = [=](const CMarker *marker_container) -> EFaceLocation
		{
			// Select the face direction based on the first element index.
			EFaceLocation iface = marker_container->GetElementFaces(0).mFace; 
			
			// Loop over each element and ensure the face is constant.
			for( auto& marker: marker_container->GetElementFaces() )
			{
				if( marker.mFace != iface )
				{
					ERROR(marker_container->GetNameMarker() + " must have the same face direction.");
				}
			}

			// Return the face direction.
			return iface;
		};

		// Deduce the faces of each marker. Note, it suffices to consider only the 
		// first element, as the entire marker must have the same face direction.
		iFaceLocation = lFaceLocation(imarker_container); 
		jFaceLocation = lFaceLocation(jmarker_container); 


		// Deduce the number of elements on both markers.
		nElem = imarker_container->GetnElem();
		
		// Check that the number of elements is not zero.
		if( nElem == 0 ) ERROR("Interface markers must not be empty.");
		
		// Ensure that both markers have the same number of elements. 
		if( nElem != jmarker_container->GetnElem() )
		{
			ERROR("Interface markers must share the same number of elements.");
		}


		// TODO: rewrite this.
		//// Initialize the elements of both pair of markers.
		//mIndexElement.reserve(mNElem);

		//for(unsigned int i=0; i<mNElem; i++)
		//{
		//	// The assumption in AS3 is that the matching face is reversed. 
		//	// This is because all zones use a clockwise convention to tag 
		//	// their boundary markers. The the matching (j-)index is:
		//	const unsigned int j = mNElem - i - 1;

		//	// Deduce the actual element indices, not their marker index.
		//	const unsigned int I = imarker_container->GetElementFaces(i).mIndex;
		//	const unsigned int J = jmarker_container->GetElementFaces(j).mIndex;

		//	mIndexElement.emplace_back(I,J);
		//}

		// TODO: adjust this.
		//// Process the markers, to ensure they coincide geometrically.
		//ProcessMatchingMarkers(config_container, 
		//		                   geometry_container, 
		//											 imarker_container, 
		//											 jmarker_container,
		//											 param_container);



	}

}

//-----------------------------------------------------------------------------------







