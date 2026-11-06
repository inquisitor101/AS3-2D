#include "geometry/face/CDirectionalFaceConnectivity.hpp"

#include <cstddef>
#include <utility>

#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "MarkerFamilies.hpp"


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
