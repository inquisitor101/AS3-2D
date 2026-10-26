#include "face_structure.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include "geometry_structure.hpp"


//-----------------------------------------------------------------------------------
// CDirectionalFaceConnectivity member functions.
//-----------------------------------------------------------------------------------

void CDirectionalFaceConnectivity::InitializeFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Function that builds the three face collections for one owner-side logical direction.
  */
{
  // Initialize the collections of faces.
  InitializeInteriorFaces(multizone_geometry_container);
  InitializeInterfaceFaces(multizone_geometry_container);
  InitializeBoundaryFaces(multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

CFlattenedFaceIndex CDirectionalFaceConnectivity::ResolveFace
(
 std::size_t i
) const
 /*
  * Function that resolves a direction-local face index into its category, batch and local face.
  */
{
  CFaceAddress face_info;
  ETypeFaceGeometry face_type;

  if( i < GetnInteriorFaces() )
  {
    face_type = ETypeFaceGeometry::INTERNAL;
    face_info = mInteriorFaces.ResolveFace(i);
  }
  else if( i < GetnInteriorFaces() + GetnBoundaryFaces() )
  {
    face_type = ETypeFaceGeometry::BOUNDARY;
    face_info = mBoundaryFaces.ResolveFace(i - GetnInteriorFaces());
  }
  else if( i < GetnFacesTotal() )
  {
    face_type = ETypeFaceGeometry::INTERFACE;
    face_info = mInterfaceFaces.ResolveFace(i - GetnInteriorFaces() - GetnBoundaryFaces());
  }
  else ERROR("Invalid direction-local face index.");

  return CFlattenedFaceIndex{ face_info.mIndexBatch, face_info.mIndexFace, face_type };
}

//-----------------------------------------------------------------------------------

void CDirectionalFaceConnectivity::InitializeInteriorFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that collects same-zone interior faces in the selected logical direction.
	*/
{
  // Temporary internal face batches, will be moved inside their collection. 
  as3vector1d<CInteriorFaceBatch> interior_batches;

  // Extract the total number of zones.
  const auto nZone = multizone_geometry_container->GetnZone();

  // Initialize the batches of interior faces, based on the number of zones.
  interior_batches.reserve( nZone );

  // Loop over each zone and initialize its internal faces, based on the class's face type.
  for( const auto& zone_geometry : multizone_geometry_container->GetSinglezoneGeometry() )
  {
    // Initialize a batch of interior faces belonging to this zone.
    interior_batches.emplace_back( mFaceDirection, zone_geometry.get() );
  }

  // Move the ownership to their respective collection.
  mInteriorFaces.Initialize( std::move(interior_batches) );
}

//-----------------------------------------------------------------------------------

void CDirectionalFaceConnectivity::InitializeInterfaceFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Function that collects internal and periodic interfaces by owner-side logical direction.
  */
{
  // Temporary interface face batches, will be moved inside their collection. 
  as3vector1d<CInterfaceFaceBatch> interface_batches;

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
        if( marker.GetiTypeFace() == mFaceDirection )
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

  // Allocate memory for the batches of interfaces.
  interface_batches.reserve( nInterfacesConsidered );

  // Loop over each internal marker family.
  for( const auto& family : internal_markers )
  {
    for( const auto& marker : family.GetMarkers() )
    {
      // Only consider interfaces whose owner face matches this class's direction.
      if( marker.GetiTypeFace() != mFaceDirection ) continue;

      // Create a batch of interface faces belonging to this zone pair.
      interface_batches.emplace_back( multizone_geometry_container, marker );
    }
  }

  // Loop over each periodic marker family.
  for( const auto& family : periodic_markers )
  {
    for( const auto& marker : family.GetMarkers() )
    {
      // Only consider interfaces whose owner face matches this class's direction.
      if( marker.GetiTypeFace() != mFaceDirection ) continue;

      // Create a batch of interface faces belonging to this periodic zone pair.
      interface_batches.emplace_back( multizone_geometry_container, marker );
    }
  }

  // Move the ownership to their respective collection.
  // Also initialize an empty collection when no interfaces are relevant.
  mInterfaceFaces.Initialize( std::move(interface_batches) );
}

//-----------------------------------------------------------------------------------

void CDirectionalFaceConnectivity::InitializeBoundaryFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Function that collects physical boundary faces in the selected logical direction.
  */
{
  // Temporary boundary face batches, will be moved inside their collection. 
  as3vector1d<CBoundaryFaceBatch> boundary_batches;

  // Extract the boundary marker families.
  const auto& boundary_markers = multizone_geometry_container->GetExternalFamilyMarkers();

   // Helper lambda to count the boundaries relevant for this class.
  auto lCountRelevantBoundaries = [&](const auto &families) -> std::size_t
  {
    std::size_t nBoundaries = 0;

    for( const auto& family : families )
    {
      for( const auto& marker : family.GetMarkers() )
      {
        if( marker.GetTypeFace() == mFaceDirection )
        {
          nBoundaries++;
        }
      }
    }

    return nBoundaries;
  };

  // Get the number of valid boundaries for this direction.
  const std::size_t nBoundariesConsidered = lCountRelevantBoundaries(boundary_markers);

  // Allocate memory for the batches of boundaries.
  boundary_batches.reserve( nBoundariesConsidered );

  // Loop over each boundary marker family.
  for( const auto& family : boundary_markers )
  {
    for( const auto& marker : family.GetMarkers() )
    {
      // Only consider boundaries whose face matches this class's direction.
      if( marker.GetTypeFace() != mFaceDirection ) continue;

      // Extract the relevant zone index.
      const std::size_t iZone = marker.mIndexZone;

      // Obtain the relevant single zone geometry.
      const auto* zone_geometry = multizone_geometry_container->GetSinglezoneGeometry(iZone);

      // Create a batch of boundary faces belonging to this face.
      boundary_batches.emplace_back( zone_geometry, marker, family.GetName() );
    }
  }

   // Move the ownership to their respective collection.
  mBoundaryFaces.Initialize( std::move(boundary_batches) );
}


//-----------------------------------------------------------------------------------
// CInteriorFaceBatch member functions.
//-----------------------------------------------------------------------------------

CInteriorFaceBatch::CInteriorFaceBatch
(
 ETypeFace                  face_direction,
 const CSinglezoneGeometry *singlezone_geometry_container
)
  : mFaceDirection( face_direction ), mIndexZone( singlezone_geometry_container->GetZoneID() )
 /*
  * Constructor that builds interior faces for one zone and logical direction.
  */
{
  InitializeFaces(singlezone_geometry_container);
}

//-----------------------------------------------------------------------------------

void CInteriorFaceBatch::InitializeFaces
(
 const CSinglezoneGeometry *singlezone_geometry_container
)
 /*
  * Function that builds same-zone interior faces and their shared neighbor offset.
  */
{        
  // Extract the relevant information for an internal face.
  const std::size_t niElem = singlezone_geometry_container->GetniElem();
  const std::size_t njElem = singlezone_geometry_container->GetnjElem();

  // Initialize the faces, based on the respective class direction.
  switch( mFaceDirection )
  {
    case(ETypeFace::IFACE):
    {
      // Deduce the number of internal i-faces in this zone.
      const std::size_t nInternalFaces = njElem * (niElem-1);

      // The offset of the right element (w.r.t. to the left) is always a +1.
      mOffsetElementPlus = 1;

			// Reserve the needed amount of memory for the internal i-faces.
			mFaces.reserve(nInternalFaces);

      // Loop over the i-faces in this zone.
      for(std::size_t j=0; j<njElem; j++)
      {
        for(std::size_t i=1; i<niElem; i++)
        {
          // Calculate the left element index, which is the "minus" element for this i-face.
          const std::size_t iElemL = j * niElem + i-1;

          // Construct the current internal i-face.
          mFaces.emplace_back( iElemL ); 
        }
      }

      // Consistency check.
      if( mFaces.size() != nInternalFaces )
      {
        ERROR("Size of the internal i-faces is wrong.");
      } 

      break;
    }

    case(ETypeFace::JFACE):
    {
      // Deduce the number of internal j-faces in this zone.
      const std::size_t nInternalFaces = niElem * (njElem-1);

      // The offset of the top element (w.r.t. to the bottom) is always a +niElem.
      mOffsetElementPlus = niElem;

			// Reserve the needed amount of memory for the internal j-faces.
			mFaces.reserve(nInternalFaces);

  		for(std::size_t j=1; j<njElem; j++)
  		{
  		  for(std::size_t i=0; i<niElem; i++)
  		  {
  		    // Deduce the flattened element indices.
  		    const std::size_t iElemB = j * niElem + i - niElem;
  		    
					mFaces.emplace_back( iElemB ); 
  		  }
  		}

      // Consistency check.
      if( mFaces.size() != nInternalFaces )
      {
        ERROR("Size of the internal j-faces is wrong.");
      } 

      break;
    }

    default: ERROR("Unknown face type encountered.");
  }
}



//-----------------------------------------------------------------------------------
// CInterfaceFaceBatch member functions.
//-----------------------------------------------------------------------------------

CInterfaceFaceBatch::CInterfaceFaceBatch
(
 const CMultizoneGeometry                     *multizone_geometry_container,
 const CInternalFamilyMarker::CInternalMarker &marker 
)
  : mIndexZoneOwner( marker.mIndexZoneI ),
    mIndexZoneNeighbor( marker.mIndexZoneJ ),
    mFaceDirectionOwner( marker.GetiTypeFace() ),
    mFaceDirectionNeighbor( marker.GetjTypeFace() ),
    mFaceLocationOwner( marker.mFaceLocationI ),
    mFaceLocationNeighbor( marker.mFaceLocationJ ),
    mIsReversed( marker.mIsReversed ),
    mTranslationVector( {C_ZERO, C_ZERO} )
 /*
  * Constructor that pairs connected zone edges and records their orientation.
  */
{
	// Deduce whether these faces are on a max (imax/jmax) or not.
	mIsMaxFaceOwner 
    = ( mFaceLocationOwner == EFaceLocation::IMAX || mFaceLocationOwner == EFaceLocation::JMAX ) ? true : false;
	mIsMaxFaceNeighbor 
    = ( mFaceLocationNeighbor == EFaceLocation::IMAX || mFaceLocationNeighbor == EFaceLocation::JMAX ) ? true : false;

  // Initiailize the interface faces belonging to these markers.
  InitializeFaces(multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

CInterfaceFaceBatch::CInterfaceFaceBatch
(
 const CMultizoneGeometry                     *multizone_geometry_container,
 const CPeriodicFamilyMarker::CPeriodicMarker &marker 
)
  : mIndexZoneOwner( marker.mIndexZoneI ),
    mIndexZoneNeighbor( marker.mIndexZoneJ ),
    mFaceDirectionOwner( marker.GetiTypeFace() ),
    mFaceDirectionNeighbor( marker.GetjTypeFace() ),
    mFaceLocationOwner( marker.mFaceLocationI ),
    mFaceLocationNeighbor( marker.mFaceLocationJ ),
    mFaceNameOwner( marker.mNameI ),
    mFaceNameNeighbor( marker.mNameJ ),
    mIsReversed( marker.mIsReversed ),
    mTranslationVector( marker.mTranslationVector )
 /*
  * Constructor that pairs periodic zone edges and records their translation.
  */
{
	// Deduce whether these faces are on a max (imax/jmax) or not.
	mIsMaxFaceOwner 
    = ( mFaceLocationOwner == EFaceLocation::IMAX || mFaceLocationOwner == EFaceLocation::JMAX ) ? true : false;
	mIsMaxFaceNeighbor 
    = ( mFaceLocationNeighbor == EFaceLocation::IMAX || mFaceLocationNeighbor == EFaceLocation::JMAX ) ? true : false;


  // Initiailize the interface faces belonging to these markers.
  InitializeFaces(multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

void CInterfaceFaceBatch::InitializeFaces
(
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Function that pairs owner and neighbor elements in the interface traversal order.
  */
{
  // Extract the local singlezone grids.
  const auto* igrid = multizone_geometry_container->GetSinglezoneGeometry( mIndexZoneOwner );
  const auto* jgrid = multizone_geometry_container->GetSinglezoneGeometry( mIndexZoneNeighbor );

  // Ensure the number of elements match.
  const std::size_t nei 
    = ( mFaceDirectionOwner == ETypeFace::IFACE ) ? igrid->GetnjElem() : igrid->GetniElem(); 
  const std::size_t nej 
    = ( mFaceDirectionNeighbor == ETypeFace::IFACE ) ? jgrid->GetnjElem() : jgrid->GetniElem();

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
  const auto ivec_ielem = igrid->ComputeSurfaceElementIndices( mFaceLocationOwner );
  const auto jvec_ielem = jgrid->ComputeSurfaceElementIndices( mFaceLocationNeighbor );

  // Consistency check.
  if( ivec_ielem.size() != nei || jvec_ielem.size() != nej )
  {
    ERROR("Number of surface elements do not match.");
  }


  // Reserve memory for the elements on this interface.
  mFaces.reserve(nbElem);

  // Loop over all the elements on this marker.
  for(std::size_t i=0; i<nbElem; i++)
  {
    // Obtain the neighbor element index, based on  whether it is reversed or not.
  	const std::size_t j = ( mIsReversed ) ? nbElem - i - 1 : i;
  
  	// Deduce the actual element indices, not their marker index.
  	const std::size_t iElem = ivec_ielem[i];
  	const std::size_t jElem = jvec_ielem[j];
  
    // Initialize the current element interface.
    mFaces.emplace_back( iElem, jElem );
  }
   
  // Ensure conformity of the markers.
  CheckConformity( multizone_geometry_container ); 
}

//-----------------------------------------------------------------------------------

void CInterfaceFaceBatch::CheckConformity
(
 const CMultizoneGeometry *multizone_geometry_container
) const
 /*
	* Function that checks paired face coordinates, orientation and periodic translation.
	*/
{
  // Abbreviations.
  const std::size_t iZone = mIndexZoneOwner;
  const std::size_t jZone = mIndexZoneNeighbor;
  
  const as3double tx = mTranslationVector[0];
  const as3double ty = mTranslationVector[1];

  // Extract the grid geometry in each of these zones.
	const auto* igrid = multizone_geometry_container->GetSinglezoneGeometry(iZone);
	const auto* jgrid = multizone_geometry_container->GetSinglezoneGeometry(jZone);

  // Get the node ndices on the ith and jth zone.
  const auto ivec_inode = igrid->ComputeSurfaceNodeIndices( mFaceLocationOwner );
  const auto jvec_inode = jgrid->ComputeSurfaceNodeIndices( mFaceLocationNeighbor );

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
  for( const auto& interface : mFaces )
  {
    // Abbreviations for the element indices.
    const std::size_t ie = interface.mIndexElementOwner;
    const std::size_t je = interface.mIndexElementNeighbor;

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


//-----------------------------------------------------------------------------------
// CBoundaryFaceBatch member functions.
//-----------------------------------------------------------------------------------

CBoundaryFaceBatch::CBoundaryFaceBatch
(
 const CSinglezoneGeometry                    *singlezone_geometry_container,
 const CExternalFamilyMarker::CExternalMarker &marker,
 const std::string                            &name
)
  : mFaceDirection( marker.GetTypeFace() ),
    mFaceLocation( marker.mFaceLocation ), 
    mIndexZone( marker.mIndexZone ),
    mName( name )
 /*
  * Constructor that builds the faces of one named physical boundary.
  */
{
  // Consistency check.
  if( mIndexZone != singlezone_geometry_container->GetZoneID() )
  {
    ERROR("Boundary has different zones!");
  }

  // Deduce whether these faces are on a max (imax/jmax) or not.
	mIsMaxFace 
    = ( mFaceLocation == EFaceLocation::IMAX || mFaceLocation == EFaceLocation::JMAX ) ? true : false;

  // Initialize the boundary faces.
  InitializeFaces(singlezone_geometry_container);
}

//-----------------------------------------------------------------------------------

void CBoundaryFaceBatch::InitializeFaces
(
 const CSinglezoneGeometry *singlezone_geometry_container
)
 /*
  * Function that enumerates the element faces along one physical boundary edge.
  */
{
  // Extract the relevant information for an external face.
  const std::size_t niElem = singlezone_geometry_container->GetniElem();
  const std::size_t njElem = singlezone_geometry_container->GetnjElem();

  // Initialize the faces, based on the respective class direction.
  switch( mFaceDirection )
  {
    case(ETypeFace::IFACE):
    {
			// Reserve the needed amount of memory for the boundary i-faces.
			mFaces.reserve(njElem);

      // Deduce the iElem index.
      const std::size_t i = ( mFaceLocation == EFaceLocation::IMIN ) ? 0 : niElem-1;

      // Safety check.
      if( mFaceLocation != EFaceLocation::IMIN && mFaceLocation != EFaceLocation::IMAX )
      {
        ERROR("Face location is inconsistent.");
      } 

      // Loop over the i-faces in this zone.
      for(std::size_t j=0; j<njElem; j++)
      {
        // Deduce the owner flattened element index of this face.
        const std::size_t ijElem = j * niElem + i;

        // Construct the current boundary i-face.
        mFaces.emplace_back( ijElem ); 
      }

      break;
    }

    case(ETypeFace::JFACE):
    {
			// Reserve the needed amount of memory for the boundary j-faces.
			mFaces.reserve(niElem);

      // Deduce the jElem index.
      const std::size_t j = ( mFaceLocation == EFaceLocation::JMIN ) ? 0 : njElem-1;

      // Safety check.
      if( mFaceLocation != EFaceLocation::JMIN && mFaceLocation != EFaceLocation::JMAX )
      {
        ERROR("Face location is inconsistent.");
      } 

      // Loop over the j-faces in this zone.
  		for(std::size_t i=0; i<niElem; i++)
  		{
  		  // Deduce the owner flattened element index of this face.
  		  const std::size_t ijElem = j * niElem + i;
  		
        // construct the current boundary j-face.
				mFaces.emplace_back( ijElem ); 
  		}

      break;
    }

    default: ERROR("Unknown face type encountered.");
  }
}







