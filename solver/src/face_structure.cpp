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

  // TODO: initialize the family of boundary faces.
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
  * only doing so for ith components whose face direction matches this class.
  * Both internal and periodic markers define interfaces.
  */
{
  // TODO: make the group initialize this, since it owns them.
  // Temporary interface families, will be moved inside their group. 
  as3vector1d<CInterfaceFacesFamily> interface_families;

  // Extract the internal and periodic marker families.
  const auto& internal_markers = multizone_geometry_container->GetInternalFamilyMarkers();
  const auto& periodic_markers = multizone_geometry_container->GetPeriodicFamilyMarkers();

   
  // Helper lambda to count the interfaces relevant for this class.
  auto lCountRelevantInterfaces = [&](const auto &families) -> std::size_t
  {
    std::size_t nInterfaces = 0;

    for( const auto& family : families )
    {
      for( const auto& marker : family.GetMarkers() )
      {
        if( marker.GetiTypeFace() == mTypeFace )
        {
          nInterfaces++;
        }
      }
    }

    return nInterfaces;
  };

  // Get the number of valid interfaces for this direction.
  const std::size_t nInterfacesConsidered = lCountRelevantInterfaces(internal_markers) 
                                          + lCountRelevantInterfaces(periodic_markers);

  // Allocate memory for the family of interfaces.
  interface_families.reserve( nInterfacesConsidered );

  // Loop over each internal marker family.
  for( const auto& family : internal_markers )
  {
    for( const auto& marker : family.GetMarkers() )
    {
      // Only consider interfaces whose ith face matches this class's direction.
      if( marker.GetiTypeFace() != mTypeFace ) continue;

      // Create a family of interface faces belonging to this zone pair.
      interface_families.emplace_back( multizone_geometry_container, marker );
    }
  }

  // Loop over each periodic marker family.
  for( const auto& family : periodic_markers )
  {
    for( const auto& marker : family.GetMarkers() )
    {
      // Only consider interfaces whose ith face matches this class's direction.
      if( marker.GetiTypeFace() != mTypeFace ) continue;

      // Create a family of interface faces belonging to this periodic zone pair.
      interface_families.emplace_back( multizone_geometry_container, marker );
    }
  }

  // Move the ownership to their respective group.
  // Also initialize an empty group when no interfaces are relevant.
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
 const CMultizoneGeometry                     *multizone_geometry_container,
 const CInternalFamilyMarker::CInternalMarker &marker 
)
  : mIndexZoneI( marker.mIndexZoneI ),
    mIndexZoneJ( marker.mIndexZoneJ ),
    mTypeFaceI( marker.GetiTypeFace() ),
    mTypeFaceJ( marker.GetjTypeFace() ),
    mFaceLocationI( marker.mFaceLocationI ),
    mFaceLocationJ( marker.mFaceLocationJ ),
    mIsReversed( marker.mIsReversed ),
    mIsPeriodic( false ),
    mTranslationVector( {C_ZERO, C_ZERO} )
 /*
  *
  */
{
	// Deduce whether these faces are on a max (imax/jmax) or not.
	mIsMaxFaceI 
    = ( mFaceLocationI == EFaceLocation::IMAX || mFaceLocationI == EFaceLocation::JMAX ) ? true : false;
	mIsMaxFaceJ 
    = ( mFaceLocationJ == EFaceLocation::IMAX || mFaceLocationJ == EFaceLocation::JMAX ) ? true : false;


  // Initiailize the interface faces belonging to these markers.
  InitializeInterfaceFaces(multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

CInterfaceFacesFamily::CInterfaceFacesFamily
(
 const CMultizoneGeometry                     *multizone_geometry_container,
 const CPeriodicFamilyMarker::CPeriodicMarker &marker 
)
  : mIndexZoneI( marker.mIndexZoneI ),
    mIndexZoneJ( marker.mIndexZoneJ ),
    mTypeFaceI( marker.GetiTypeFace() ),
    mTypeFaceJ( marker.GetjTypeFace() ),
    mFaceLocationI( marker.mFaceLocationI ),
    mFaceLocationJ( marker.mFaceLocationJ ),
    mIsReversed( marker.mIsReversed ),
    mFaceNameI( marker.mNameI ),
    mFaceNameJ( marker.mNameJ ),
    mIsPeriodic( true ),
    mTranslationVector( marker.mTranslationVector )
 /*
  *
  */
{
	// Deduce whether these faces are on a max (imax/jmax) or not.
	mIsMaxFaceI 
    = ( mFaceLocationI == EFaceLocation::IMAX || mFaceLocationI == EFaceLocation::JMAX ) ? true : false;
	mIsMaxFaceJ 
    = ( mFaceLocationJ == EFaceLocation::IMAX || mFaceLocationJ == EFaceLocation::JMAX ) ? true : false;


  // Initiailize the interface faces belonging to these markers.
  InitializeInterfaceFaces(multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

void CInterfaceFacesFamily::InitializeInterfaceFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  *
  */
{
  // Extract the local singlezone grids.
  const auto* igrid = multizone_geometry_container->GetSinglezoneGeometry( mIndexZoneI );
  const auto* jgrid = multizone_geometry_container->GetSinglezoneGeometry( mIndexZoneJ );

  // Ensure the number of elements match.
  const std::size_t nei 
    = ( mTypeFaceI == ETypeFace::IFACE ) ? igrid->GetnjElem() : igrid->GetniElem(); 
  const std::size_t nej 
    = ( mTypeFaceJ == ETypeFace::IFACE ) ? jgrid->GetnjElem() : jgrid->GetniElem();

  // Consistency check.
  if( nei != nej )
  {
    ERROR("Number of elements on this interface do not match.");
  }

  // Select one as the number of elements on this boundary.
  const std::size_t nbElem = nei; 
  
  // Check that the number of elements is not zero.
  if( nbElem == 0 ) ERROR("Interface markers must not be empty.");
 

  // Get the element indices on this face.
  const auto ivec_ielem = igrid->ComputeSurfaceElementIndices( mFaceLocationI );
  const auto jvec_ielem = jgrid->ComputeSurfaceElementIndices( mFaceLocationJ );

  // Consistency check.
  if( ivec_ielem.size() != nei || jvec_ielem.size() != nej )
  {
    ERROR("Number of surface elements do not match.");
  }


  // Reserve memory for the elements on this interface.
  mInterfaceFaces.reserve(nbElem);

  // Loop over all the elements on this marker.
  for(std::size_t i=0; i<nbElem; i++)
  {
    // Obtain the neighbor element index, based on  whether it is reversed or not.
  	const std::size_t j = ( mIsReversed ) ? nbElem - i - 1 : i;
  
  	// Deduce the actual element indices, not their marker index.
  	const std::size_t iElem = ivec_ielem[i];
  	const std::size_t jElem = jvec_ielem[j];
  
    // Initialize the current element interface.
    mInterfaceFaces.emplace_back( iElem, jElem );
  }
   
  // Ensure conformity of the markers.
  CheckConformityMarkers( multizone_geometry_container ); 
}

//-----------------------------------------------------------------------------------

void CInterfaceFacesFamily::CheckConformityMarkers
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that processes each pair of markers, such that their common face coincides.
	*/
{
  // Abbreviations.
  const std::size_t iZone = mIndexZoneI;
  const std::size_t jZone = mIndexZoneJ;
  
  const as3double tx = mTranslationVector[0];
  const as3double ty = mTranslationVector[1];

  // Extract the grid geometry in each of these zones.
	const auto* igrid = multizone_geometry_container->GetSinglezoneGeometry(iZone);
	const auto* jgrid = multizone_geometry_container->GetSinglezoneGeometry(jZone);

  // Get the node ndices on the ith and jth zone.
  const auto ivec_inode = igrid->ComputeSurfaceNodeIndices( mFaceLocationI );
  const auto jvec_inode = jgrid->ComputeSurfaceNodeIndices( mFaceLocationJ );

  // Report an error, in case mismatch in number of nodes.
  if( ivec_inode.size() != jvec_inode.size() )
  {
    ERROR("Polynomial order must be the same in the two zones of a periodic boundary.");
  }

  // We only check the endpoints on an element for conformity.
  const std::size_t in0 = ivec_inode.front();
  const std::size_t in1 = ivec_inode.back();

  // Swap the j endpoints when checking the reversed orientation.
  const std::size_t jn0 = mIsReversed ? jvec_inode.back()  : jvec_inode.front();
  const std::size_t jn1 = mIsReversed ? jvec_inode.front() : jvec_inode.back();

  // Relative tolerance value.
  const as3double rtol = static_cast<as3double>( 1.0e-8 );

  // Absolute tolerance value, in grid coordinate units.
  const as3double atol = static_cast<as3double>( 1.0e-10 );

  // Compare a translated i coordinate with its matching j coordinate.
  auto lCoordinatesMatch = [rtol, atol](as3double xi, as3double xj)
  {
    const as3double difference = std::abs( xi - xj );
    const as3double tolerance = atol + rtol * std::max( std::abs(xi), std::abs(xj) );

    return difference <= tolerance;
  };


  // Loop over every element on this interface.
  for( const auto& interface : mInterfaceFaces )
  {
    // Abbreviations for the element indices.
    const std::size_t ie = interface.mIndexElementI;
    const std::size_t je = interface.mIndexElementJ;

    // Extract the entire element coordinates.
    const auto& icoor = igrid->GetElementGeometry( ie )->GetCoordSolDOFs();
    const auto& jcoor = jgrid->GetElementGeometry( je )->GetCoordSolDOFs();

    // Extract the endpoint coordinates for the current element edge, apply the translation.
    const as3double xi0 = icoor( 0, in0 ) + tx, xi1 = icoor( 0, in1 ) + tx;
    const as3double yi0 = icoor( 1, in0 ) + ty, yi1 = icoor( 1, in1 ) + ty;

    // Extract coordinates for the matching marker.
    // For the reversed case, jn0 and jn1 already select the opposite endpoints.
    const as3double xj0 = jcoor( 0, jn0 ), xj1 = jcoor( 0, jn1 );
    const as3double yj0 = jcoor( 1, jn0 ), yj1 = jcoor( 1, jn1 );

    // If either endpoint does not match, report error.
    if( !lCoordinatesMatch(xi0, xj0) || !lCoordinatesMatch(yi0, yj0) ||
        !lCoordinatesMatch(xi1, xj1) || !lCoordinatesMatch(yi1, yj1) )
    {
      ERROR( "Interface does not match with "
             " iZone: " + std::to_string(iZone) + ","
             " jZone: " + std::to_string(jZone) ); 
    }
  }
}





