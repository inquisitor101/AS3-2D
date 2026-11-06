#include "geometry/face/CInterfaceFaceBatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>

#include "geometry/CElementGeometry.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"


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
