#include "geometry/face/CInteriorFaceBatch.hpp"

#include <cstddef>

#include "geometry/CSinglezoneGeometry.hpp"


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
