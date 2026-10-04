#include "face_structure.hpp"
#include "geometry_structure.hpp"


//-----------------------------------------------------------------------------------
// CMultizoneFaceGeometry member functions.
//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::InitializeFaces
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  *
  */
{
  // Initialize the family of faces.
  InitializeInternalFaces(multizone_geometry_container);
  InitializeInterfaceFaces(config_container, multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::InitializeInternalFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	*
	*/
{
  // TODO: make the group initialize this, since it owns them.
  // Temporary internal families, will be moved inside their group. 
  as3vector1d<CInternalFacesFamily> internal_families;

  // Extract the total number of zones.
  const auto nZone = multizone_geometry_container->GetnZone();

  // Initialize the family of internal faces, based on the number of zones.
  internal_families.reserve( nZone );

  // Loop over each zone and initialize its internal faces, based on the class's face type.
  for( const auto& zone_geometry : multizone_geometry_container->GetSinglezoneGeometry() )
  {
    // Initialize a family of internal faces belonging to this zone.
    internal_families.emplace_back( mTypeFace, zone_geometry.get() );
  }

  // Move the ownership to their respective group.
  mInternalFacesGroup.Initialize( std::move(internal_families) );
}

//-----------------------------------------------------------------------------------

void CMultizoneFaceGeometry::InitializeInterfaceFaces
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that initializes interfaces across all zones with the assumption of
  * only doing so for ith components who have faces in the i-direction.
	*/
{
  // TODO: make the group initialize this, since it owns them.
  // Temporary interface families, will be moved inside their group. 
  as3vector1d<CInterfaceFacesFamily> interface_families;

  // Helper lambda to find the matching marker.
  auto lGetMatchingMarker = [](const CMultizoneGeometry *multizone_geometry_container, 
                               const std::string &marker_name) -> const CMarker*
  {
    for(auto& zone : multizone_geometry_container->GetSinglezoneGeometry())
    {
      for(auto& marker : zone->GetMarker())
      {
        if(marker->GetNameMarker() == marker_name) return marker.get();
      }
    }
  
    ERROR("Could not find the interface marker.");
  };


  // Helper lambda to determinwe which interface is relevant for this class.
  auto lGetRelevantIndexInterfaces = [&](const CConfig *config_container, 
                                         const CMultizoneGeometry *multizone_geometry_container) -> as3vector1d<size_t>
  {
    // Variable for book-keeping indices of relevant interfaces.
    as3vector1d<size_t> interfaces_considered;
  
    // Extract all the interfaces specified by the user.
    auto& param_interface_total = config_container->GetInterfaceParamMarker();
 
    // If there are no interfaces, we leave.
    if( param_interface_total.empty() ) return interfaces_considered;
 
    // Loop over all interfaces and check which to include (i.e. are relevant).
    for(size_t i=0; i<param_interface_total.size(); i++)
    {
      // Extract the current interface.
      const auto& param_interface = param_interface_total[i];
  
      // Extract its marker container.
      const CMarker* imarker_container = lGetMatchingMarker(multizone_geometry_container, param_interface->mName);
  
      // If this face's ith face type matches the one assigned for this class, add it.
      if( imarker_container->GetTypeFace() == mTypeFace )
      {
        interfaces_considered.emplace_back(i);
      }
    }
  
    return interfaces_considered;
  };




	// Extract the interface boundaries, specified by the user.
	auto& param_interface_total = config_container->GetInterfaceParamMarker();

  // Get the valid interfaces, based on the ith face and the direction in this class.
	const as3vector1d<size_t> interfaces_considered = lGetRelevantIndexInterfaces(config_container, 
                                                                                multizone_geometry_container);
 
	// Get the number of valid interfaces for this direction.
  const size_t nInterfacesConsidered = interfaces_considered.size();

	// If there are no valid interfaces, return without any initialization.
	if( nInterfacesConsidered == 0 ) return;

	// Allocate memory for the family of interfaces.
	interface_families.reserve( nInterfacesConsidered );

  // Loop over each interface group.
  for( const size_t i : interfaces_considered )
  {
    // Extract the relevant marker interface.
    const auto& param_interface = param_interface_total[i];

    // Get a pointer to the two markers forming this interface.
    const CMarker *imarker_container = lGetMatchingMarker(multizone_geometry_container, param_interface->mName);
    const CMarker *jmarker_container = lGetMatchingMarker(multizone_geometry_container, param_interface->mNameMatching);

    // Create a family of interface faces belonging to these markers.
    interface_families.emplace_back( multizone_geometry_container, 
                                     imarker_container, 
                                     jmarker_container, 
                                     param_interface.get() );
  }
  
  // Move the ownership to their respective group.
  mInterfaceFacesGroup.Initialize( std::move(interface_families) );
}


//-----------------------------------------------------------------------------------
// CInternalFacesFamily member functions.
//-----------------------------------------------------------------------------------

CInternalFacesFamily::CInternalFacesFamily
(
 ETypeFace                  face_type,
 const CSinglezoneGeometry *singlezone_geometry_container
)
  : mTypeFace( face_type ), mIndexZone( singlezone_geometry_container->GetZoneID() )
 /*
  *
  */
{
  InitializeInternalFaces(singlezone_geometry_container);
}

//-----------------------------------------------------------------------------------

void CInternalFacesFamily::InitializeInternalFaces
(
 const CSinglezoneGeometry *singlezone_geometry_container
)
 /*
  *
  */
{        
  // Extract the relevant information for an internal face.
  const size_t niElem = singlezone_geometry_container->GetniElem();
  const size_t njElem = singlezone_geometry_container->GetnjElem();

  // Initialize the faces, based on the respective class direction.
  switch( mTypeFace )
  {
    case(ETypeFace::IFACE):
    {
      // Deduce the number of internal i-faces in this zone.
      const size_t nInternalFaces = njElem * (niElem-1);

      // The offset of the right element (w.r.t. to the left) is always a +1.
      mOffsetElementPlus = 1;

			// Reserve the needed amount of memory for the internal i-faces.
			mInternalFaces.reserve(nInternalFaces);

      // Loop over the i-faces in this zone.
      for(size_t j=0; j<njElem; j++)
      {
        for(size_t i=1; i<niElem; i++)
        {
          // Calculate the left element index, which is the "minus" element for this i-face.
          const size_t iElemL = j * niElem + i-1;

          // Construct the current internal i-face.
          mInternalFaces.emplace_back( iElemL ); 
        }
      }

      // Consistency check.
      if( mInternalFaces.size() != nInternalFaces )
      {
        ERROR("Size of the internal i-faces is wrong.");
      } 

      break;
    }

    case(ETypeFace::JFACE):
    {
      // Deduce the number of internal j-faces in this zone.
      const size_t nInternalFaces = niElem * (njElem-1);

      // The offset of the top element (w.r.t. to the bottom) is always a +niElem.
      mOffsetElementPlus = niElem;

			// Reserve the needed amount of memory for the internal j-faces.
			mInternalFaces.reserve(nInternalFaces);

  		for(size_t j=1; j<njElem; j++)
  		{
  		  for(size_t i=0; i<niElem; i++)
  		  {
  		    // Deduce the flattened element indices.
  		    const size_t iElemB = j * niElem + i - niElem;
  		    
					mInternalFaces.emplace_back( iElemB ); 
  		  }
  		}

      // Consistency check.
      if( mInternalFaces.size() != nInternalFaces )
      {
        ERROR("Size of the internal j-faces is wrong.");
      } 

      break;
    }

    default: ERROR("Unknown face type encountered.");
  }
}



//-----------------------------------------------------------------------------------
// CInterfaceFacesFamily member functions.
//-----------------------------------------------------------------------------------

CInterfaceFacesFamily::CInterfaceFacesFamily
(
 const CMultizoneGeometry *multizone_geometry_container,
 const CMarker            *imarker_container,
 const CMarker            *jmarker_container,
 CInterfaceParamMarker    *param_interface
)
 /*
  *
  */
{
  // Initiailize the interface faces belonging to these markers.
  InitializeInterfaceFaces(multizone_geometry_container,
                           imarker_container,
                           jmarker_container,
                           param_interface);
}

//-----------------------------------------------------------------------------------

void CInterfaceFacesFamily::InitializeInterfaceFaces
(
 const CMultizoneGeometry *multizone_geometry_container,
 const CMarker            *imarker_container,
 const CMarker            *jmarker_container,
 CInterfaceParamMarker    *param_interface
)
 /*
  *
  */
{
  // Extract the zone ID of these markers.
  mIndexZoneI = imarker_container->GetZoneID();
  mIndexZoneJ = jmarker_container->GetZoneID();

  // Extract the face location of these markers.
  mFaceLocationI = imarker_container->GetFaceLocation(); 
  mFaceLocationJ = jmarker_container->GetFaceLocation(); 
 
  // Extract the face names of these markers.
  mFaceNameI = imarker_container->GetNameMarker();
  mFaceNameJ = jmarker_container->GetNameMarker();

	// Deduce whether these faces are on a max (imax/jmax) or not.
	mIsMaxFaceI = ( mFaceLocationI == EFaceLocation::IMAX || mFaceLocationI == EFaceLocation::JMAX ) ? true : false;
	mIsMaxFaceJ = ( mFaceLocationJ == EFaceLocation::IMAX || mFaceLocationJ == EFaceLocation::JMAX ) ? true : false;

  // Deduce the number of elements on both markers.
  const size_t nElem = imarker_container->GetnElem();
  
  // Check that the number of elements is not zero.
  if( nElem == 0 ) ERROR("Interface markers must not be empty.");
  
  // Ensure that both markers have the same number of elements. 
  if( nElem != jmarker_container->GetnElem() )
  {
  	ERROR("Interface markers must share the same number of elements.");
  }
  
  // Reserve memory for the elements on this interface.
  mInterfaceFaces.reserve(nElem);

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
    mInterfaceFaces.emplace_back( iElem, jElem );
  }
   
  // Ensure conformity of the markers.
  CheckConformityMarkers( multizone_geometry_container, 
  											  imarker_container, 
  											  jmarker_container,
  											  param_interface ); 
}

//-----------------------------------------------------------------------------------

void CInterfaceFacesFamily::CheckConformityMarkers
(
 const CMultizoneGeometry *multizone_geometry_container,
 const CMarker            *imarker_container,
 const CMarker            *jmarker_container,
 CInterfaceParamMarker    *param_interface
)
 /*
	* Function that processes each pair of markers, such that their common face coincides.
	*/
{
  // Extract the relevant information in this marker.
  const unsigned short iZone = mIndexZoneI;
  const unsigned short jZone = mIndexZoneJ;
  const auto&  iName         = mFaceNameI;
  const auto&  jName         = mFaceNameJ;
  const size_t nFace         = GetnFace();
	
  // Extract the grid geometry in each of these zones.
	auto* igrid = multizone_geometry_container->GetSinglezoneGeometry(iZone);
	auto* jgrid = multizone_geometry_container->GetSinglezoneGeometry(jZone);

	// Extract the properties of each marker region (element indices and faces).
	auto& imarker = imarker_container->GetElementFaces(); 
	auto& jmarker = jmarker_container->GetElementFaces();

	// Ensure the number of indices in each marker matches.
	if( (imarker.size() != jmarker.size()) || (imarker.size() != nFace) ) 
	{
		ERROR("Interface markers have different number of elements.");
	}

	// Relative tolerance value.
	const as3double tol = static_cast<as3double>( 1.0e-8 );

	// Loop over each pair of elements on this face and check that their faces coincide.
	for(size_t i=0; i<nFace; i++)
	{
		// The assumption in AS3 is that the matching face is reversed. 
		// This is because all zones use a clockwise convention to tag 
		// their boundary markers. The the matching (j-)index is:
		const size_t j = nFace - i - 1;

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
    const size_t iElem = mInterfaceFaces[i].mIndexElementI;
    const size_t jElem = mInterfaceFaces[i].mIndexElementJ;

    // Additional consistency check.
    if ( (iElem != imarker[i].mIndex) or (jElem != jmarker[j].mIndex) )
    {
      ERROR("Interface element indices do not coincide on a shared face.");
    }
  }
}





