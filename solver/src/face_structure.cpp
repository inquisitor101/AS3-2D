#include "face_structure.hpp"
#include "geometry_structure.hpp"


//-----------------------------------------------------------------------------------
// CZoneGeometry member functions.
//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::InitializeInternalFacesIDir
(
 const CGeometry *geometry_container
)
 /*
	*
	*/
{
  // Consistency check.
  if( mDirection != ETypeDirection::IDIR ) ERROR("Wrong function called.");

	// Determine total number of internal faces in the i-direction.
	size_t nInternalFacesIDir = 0;
	for( const auto& zone : geometry_container->GetZoneGeometry() )
	{
		const size_t niElem = zone->GetniElem();
    const size_t njElem = zone->GetnjElem();
	
		nInternalFacesIDir += njElem * (niElem-1);
	}
	
	// Reserve the needed amount of memory for the internal faces.
	mInternalFaces.reserve(nInternalFacesIDir);

	// Loop over the internal faces in all the zones and instantiate them.
  for( const auto& zone : geometry_container->GetZoneGeometry() )
  {
    const unsigned short iZone = zone->GetZoneID();
    const size_t niElem = zone->GetniElem();
    const size_t njElem = zone->GetnjElem();

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

void CMultizoneFaceGeometry::InitializeInternalFacesJDir
(
 const CGeometry *geometry_container
)
 /*
	*
	*/
{
  // Consistency check.
  if( mDirection != ETypeDirection::JDIR ) ERROR("Wrong function called.");

	// Determine total number of internal faces in the j-direction.
	size_t nInternalFacesJDir = 0;
	for( const auto& zone : geometry_container->GetZoneGeometry() )
	{
		const size_t niElem = zone->GetniElem();
    const size_t njElem = zone->GetnjElem();
	
		nInternalFacesJDir += niElem * (njElem-1);
	}
	
	// Reserve the needed amount of memory for the internal faces.
	mInternalFaces.reserve(nInternalFacesJDir);

	// Loop over the internal faces in all the zones and instantiate them.
  for( const auto& zone : geometry_container->GetZoneGeometry() )
  {
    const unsigned short iZone = zone->GetZoneID();
    const size_t niElem = zone->GetniElem();
    const size_t njElem = zone->GetnjElem();

    for(size_t j=1; j<njElem; j++)
    {
      for(size_t i=0; i<niElem; i++)
      {
        // Deduce the flattened element indices.
        const size_t ijElemT = j * niElem + i;
        const size_t ijElemB = ijElemT - niElem;
        
				mInternalFaces.emplace_back( ijElemB, ijElemT, iZone ); 
      }
    }
	}
}

//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::InitializeInterfaceFacesIDir
(
 const CConfig   *config_container,
 const CGeometry *geometry_container
)
 /*
	* Function that initializes interfaces across all zones with the assumption of
  * only doing so for ith components who have faces in the i-direction.
	*/
{
  // Consistency check.
  if( mDirection != ETypeDirection::IDIR ) ERROR("Wrong function called.");

	// Extract the interface boundaries, specified by the user.
	auto& param_interface_total = config_container->GetInterfaceParamMarker();

  // Get the number of valid interfaces for this direction.
  //const size_t nInterfacesConsidered = GetnInterfacesAlongDirection(config_container, geometry_container);
  const as3vector1d<size_t> interfaces_considered = GetIndexInterfacesAlongDirection(config_container, geometry_container);
  
  const size_t nInterfacesConsidered = interfaces_considered.size();

	// If there are no valid interfaces, return without any initialization.
	if( nInterfacesConsidered == 0 ) return;

	// Allocate memory for the family of interfaces.
	mInterfaceGroups.reserve( nInterfacesConsidered );

  // Loop over each interface group.
  for( const size_t i : interfaces_considered )
  {
    const auto& param_interface = param_interface_total[i];

    // Get a pointer to the two markers forming this interface.
    const CMarker *imarker_container = GetMatchingMarker(geometry_container, param_interface->mName);
    const CMarker *jmarker_container = GetMatchingMarker(geometry_container, param_interface->mNameMatching);
    
    // Extract the zone ID of these markers.
    const unsigned short iZone = imarker_container->GetZoneID();
    const unsigned short jZone = jmarker_container->GetZoneID();

    // Deduce the number of elements on both markers.
    const size_t nElem = imarker_container->GetnElem();
    
    // Check that the number of elements is not zero.
    if( nElem == 0 ) ERROR("Interface markers must not be empty.");
    
    // Ensure that both markers have the same number of elements. 
    if( nElem != jmarker_container->GetnElem() )
    {
    	ERROR("Interface markers must share the same number of elements.");
    }
 
    // Deduce the faces of each marker. Note, it suffices to consider only the 
    // first element, as the entire marker must have the same face direction.
    const EFaceLocation iFaceLocation = GetMarkerFaceLocation(imarker_container); 
    const EFaceLocation jFaceLocation = GetMarkerFaceLocation(jmarker_container); 

    // If the face location of the ith element faces are not in the i-direction, skip this.
    if( GetDirectionFromFaceLocation(iFaceLocation) != mDirection ) continue;

    // First face belonging to this interface group.
    const size_t interface_begin = mInterfaceFaces.size();

    // Loop over all the elements on this marker.
    for(size_t i=0; i<nElem; i++)
    {
    	// The assumption in AS3 is that the matching face is reversed. 
    	// This is because all zones use a clockwise convention to tag 
    	// their boundary markers. The the matching (j-)index is:
    	const size_t j = nElem - i - 1;
    
    	// Deduce the actual element indices, not their marker index.
    	const size_t iElem = imarker_container->GetElementFaces(i).mIndex;
    	const size_t jElem = jmarker_container->GetElementFaces(j).mIndex;
   
      // Initialize the current element interface.
      mInterfaceFaces.emplace_back( iElem, jElem, iZone, jZone, iFaceLocation, jFaceLocation );
    }
    
    // One group describing the entire contiguous range.
    mInterfaceGroups.emplace_back( iZone, jZone, 
                                   iFaceLocation, jFaceLocation, 
                                   imarker_container->GetNameMarker(), jmarker_container->GetNameMarker(),
                                   interface_begin, mInterfaceFaces.size() );

    // Ensure conformity of the markers. NOTE, carefull with .back(),
    // if the size changes before this call (however, this should never happen here).
    CheckConformityMarkers(config_container, 
	  		                   geometry_container, 
	  											 imarker_container, 
	  											 jmarker_container,
	  											 param_interface.get(),
                           &mInterfaceGroups.back() ); 
  }
}

//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::CheckConformityMarkers
(
 const CConfig         *config_container,
 const CGeometry       *geometry_container,
 const CMarker         *imarker_container,
 const CMarker         *jmarker_container,
 CInterfaceParamMarker *param_interface,
 CInterfaceGroup       *interface_group
)
 /*
	* Function that processes each pair of markers, such that their common face coincides.
	*/
{
  // Extract the relevant information in this marker.
  const unsigned short iZone = interface_group->mIndexZoneM;
  const unsigned short jZone = interface_group->mIndexZoneP;
  const size_t nElem         = interface_group->GetnElem();
  const auto&  iName         = interface_group->mFaceNameM;
  const auto&  jName         = interface_group->mFaceNameP;

	// Extract the grid geometry in each of these zones.
	auto* igrid = geometry_container->GetZoneGeometry(iZone);
	auto* jgrid = geometry_container->GetZoneGeometry(jZone);

	// Extract the properties of each marker region (element indices and faces).
	auto& imarker = imarker_container->GetElementFaces(); 
	auto& jmarker = jmarker_container->GetElementFaces();

	// Ensure the number of indices in each marker matches.
	if( (imarker.size() != jmarker.size()) || (imarker.size() != nElem) ) 
	{
		ERROR("Interface markers have different number of elements.");
	}

	// Relative tolerance value.
	const as3double tol = static_cast<as3double>( 1.0e-8 );

	// Loop over each pair of elements and check that their faces coincide.
	for(size_t i=0; i<nElem; i++)
	{
		// The assumption in AS3 is that the matching face is reversed. 
		// This is because all zones use a clockwise convention to tag 
		// their boundary markers. The the matching (j-)index is:
		const size_t j = nElem - i - 1;

		// Extract the nodal indices of the first element on this marker.
		auto& icoor = igrid->GetElementGeometry( imarker[i].mIndex )->GetCoordSolDOFs(); 
		auto& jcoor = jgrid->GetElementGeometry( jmarker[j].mIndex )->GetCoordSolDOFs(); 
		
		// Extract the nodal indices from the face type.
		auto& inode = igrid->GetFaceNodalIndices( imarker[i].mFace );
		auto& jnode = jgrid->GetFaceNodalIndices( jmarker[j].mFace ); 

		// NOTE 
		// For now, force the polynomial orders to be the same. Otherwise, we 
		// have to come up with an integration rule that is based on the higher
		// order polynomial, to ensure local conservation.
		if( inode.size() != jnode.size() )
		{
			ERROR("Polynomial order must be identical (for now) in zones: "
					  + std::to_string(iZone) + ", " + std::to_string(jZone));
		}

		// Check if the periodic faces coincide, after translation.
		for(size_t k=0; k<inode.size(); k++)
		{
			// Extract coordinates for the current  marker.
			const as3double ix = icoor(0, inode[k]);
			const as3double iy = icoor(1, inode[k]);
			// Extract coordinates for the matching marker.
			const as3double jx = jcoor(0, jnode[k]);
			const as3double jy = jcoor(1, jnode[k]);
			
			// Compute the difference between them in absolute.
			const as3double dx = std::abs( ix - jx + param_interface->mVectorTrans[0] );
			const as3double dy = std::abs( iy - jy + param_interface->mVectorTrans[1] );

			// Relative error, based on the values of the coordinates.
			const as3double xtol = tol*std::max( std::abs(ix), std::abs(jx) );
			const as3double ytol = tol*std::max( std::abs(iy), std::abs(jy) );

			// If the boundaries do not coincide, abort with an error.
			if( (dx > xtol) || (dy > ytol) ) 
			{
				ERROR("Interface boundaries: " + iName + ", " + jName + " do not match.");
			}
		}

    // Extract the current element index of the ith side.
    const size_t iElem = mInterfaceFaces[interface_group->GetFaceIndex(i)].mIndexElementM;
    const size_t jElem = mInterfaceFaces[interface_group->GetFaceIndex(i)].mIndexElementP;

    // Additional consistency check.
    if ( (iElem != imarker[i].mIndex) or (jElem != jmarker[j].mIndex) )
    {
      ERROR("Interface element indices do not coincide on a shared face.");
    }
  }
}

//-----------------------------------------------------------------------------------

const CMarker *CMultizoneFaceGeometry::GetMatchingMarker
(
 const CGeometry   *geometry_container,
 const std::string &marker_name
)
 /*
  *
  */
{
  for( auto& zone: geometry_container->GetZoneGeometry() )
  {
  	for( auto& marker: zone->GetMarker() )
  	{
  		if( marker->GetNameMarker() == marker_name )
  		{
  			return marker.get();
  		}
  	}
  }
  
  // The program should not reach here, otherwise we have an error.
  ERROR("Could not find the interface marker.");
  
  // To avoid compiler problems, return something.
  return nullptr;
}

//-----------------------------------------------------------------------------------

const EFaceLocation CMultizoneFaceGeometry::GetMarkerFaceLocation
(
 const CMarker *marker_container
)
 /*
  *
  */
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
}

//-----------------------------------------------------------------------------------

as3vector1d<size_t> CMultizoneFaceGeometry::GetIndexInterfacesAlongDirection
(
 const CConfig   *config_container,
 const CGeometry *geometry_container
)
 /*
  *
  */
{
  // Variable for book-keeping indices.
	as3vector1d<size_t> interfaces_considered;
  
  // Extract the interface boundaries, specified by the user.
	auto& param_interface_total = config_container->GetInterfaceParamMarker();

	// If there are no interfaces specified, return without any initialization.
	if( param_interface_total.empty() ) return interfaces_considered;

  // Loop over all possible interfaces and check how many are in the direction of this class.
  for( size_t i=0; i<param_interface_total.size(); i++)
  {
    const auto& param_interface = param_interface_total[i];

    const CMarker* imarker_container = GetMatchingMarker(geometry_container, param_interface->mName);
    
    const EFaceLocation iFaceLocation = GetMarkerFaceLocation(imarker_container);
    
    if( GetDirectionFromFaceLocation(iFaceLocation) == mDirection )
    {
      interfaces_considered.emplace_back( i );
    }
  }
  
  return interfaces_considered;
}


