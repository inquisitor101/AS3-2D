#include "geometry/face/CBoundaryFaceBatch.hpp"

#include <cstddef>
#include <string>

#include "geometry/CSinglezoneGeometry.hpp"


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
