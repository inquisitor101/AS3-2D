#include "geometry/CSinglezoneGeometry.hpp"

#include <cstddef>
#include <memory>

#include "geometry/CElementGeometry.hpp"


//-----------------------------------------------------------------------------------
// CSinglezoneGeometry member functions.
//-----------------------------------------------------------------------------------


CSinglezoneGeometry::CSinglezoneGeometry
(
  const std::string      &gridfile,
  std::size_t            iZone,
  std::size_t            nPoly,
  std::size_t            niElem,
  std::size_t            njElem,
  ETypeNodalDistribution nodalDistribution,
  bool                   isAffine
)
  : mZoneID(iZone),
    mGridFile(gridfile),
    mNPolyGrid(nPoly),
    mNiElem(niElem),
    mNjElem(njElem),
    mIsAffine(isAffine),
    mNodalDistribution(nodalDistribution)
 /*
  * Constructor that records zone properties and allocates element coordinates.
  */
{
  if( mNiElem == 0 || mNjElem == 0 )
  {
    ERROR("Invalid element dimensions: " + mGridFile);
  }

  // Geometry face operations currently require LGL or equidistant solution nodes.
  if( mNodalDistribution != ETypeNodalDistribution::LGL &&
      mNodalDistribution != ETypeNodalDistribution::EQD )
  {
    ERROR("Unsupported geometry nodal distribution: " + mGridFile);
  }

  const std::size_t nNode1D = std::size_t(mNPolyGrid) + 1;
  const std::size_t nNode2D = nNode1D * nNode1D;
  const std::size_t nElem = std::size_t(mNiElem) * mNjElem;

  // Generate the elements.
  mElementGeometry.reserve(nElem);
  for(std::size_t iElem=0; iElem<nElem; iElem++)
  {
    mElementGeometry.push_back( std::make_unique<CElementGeometry>(nNode2D) );
  }
}

//-----------------------------------------------------------------------------------

CSinglezoneGeometry::~CSinglezoneGeometry
(
 void
)
 /*
	* Destructor for CSinglezoneGeometry.
	*/
{

}

//-----------------------------------------------------------------------------------

as3vector1d<std::size_t> CSinglezoneGeometry::ComputeSurfaceElementIndices
(
 EFaceLocation face_location
) const
 /*
  * Function that returns zone-local element indices in traversal order along an edge.
  */
{
  // Prevent unsigned underflow when computing niElem-1 or njElem-1.
  if( mNiElem == 0 || mNjElem == 0 )
  {
    ERROR("Cannot compute surface indices for an empty zone.");
  }

  // Initialize a structure for common values.
  struct
  {
    std::size_t mNbElem;
    std::size_t mStride;
    std::size_t mIStart;
  } common{};

  // The below returns the element indices on a surface, based on the local 
  // orientation and local indexing (hence the combination of cw and ccw).
  // i.e.,
  //        for(jElem<njElem)
  //          for(iElem<niElem) index = jElem * niElem + iElem
  // for, 
  //  *) IMIN: iElem = 0
  //  *) IMAX: iElem = niElem-1
  //  *) JMIN: jElem = 0
  //  *) JMAX: jElem = njElem-1
  switch( face_location )
  {
    case( EFaceLocation::JMIN ): { common = {mNiElem,       1,                     0}; break; }
    case( EFaceLocation::JMAX ): { common = {mNiElem,       1, (mNjElem-1) * mNiElem}; break; }
    case( EFaceLocation::IMIN ): { common = {mNjElem, mNiElem,                     0}; break; }
    case( EFaceLocation::IMAX ): { common = {mNjElem, mNiElem,             mNiElem-1}; break; }
    default: ERROR("Unknown face location input.");
  }

  // Populate the element indices.
  as3vector1d<std::size_t> indices( common.mNbElem );
  for(std::size_t k=0; k<common.mNbElem; k++)
  {
    indices[k] = k * common.mStride + common.mIStart;
  }

  return indices;
} 

//-----------------------------------------------------------------------------------

as3vector1d<std::size_t> CSinglezoneGeometry::ComputeSurfaceNodeIndices
(
 EFaceLocation face_location
) const
 /*
  * Function that returns reference-node indices in traversal order along an edge.
  */
{
  // Number of nodes in each local direction of the element.
  const std::size_t nNode1D = mNPolyGrid + 1;

  // Initialize a structure for common values.
  struct
  {
    std::size_t mStride;
    std::size_t mIStart;
  } common{};

  // The below returns the node indices on a surface, based on the local
  // orientation and local indexing (hence the combination of cw and ccw).
  // i.e.,
  //        for(jNode<nNode1D)
  //          for(iNode<nNode1D) index = jNode * nNode1D + iNode
  // for,
  //  *) IMIN: iNode = 0
  //  *) IMAX: iNode = nNode1D-1
  //  *) JMIN: jNode = 0
  //  *) JMAX: jNode = nNode1D-1
  switch( face_location )
  {
    case( EFaceLocation::JMIN ): { common = {      1,                     0}; break; }
    case( EFaceLocation::JMAX ): { common = {      1, (nNode1D-1) * nNode1D}; break; }
    case( EFaceLocation::IMIN ): { common = {nNode1D,                     0}; break; }
    case( EFaceLocation::IMAX ): { common = {nNode1D,             nNode1D-1}; break; }
    default: ERROR("Unknown face location input.");
  }

  // Populate the node indices.
  as3vector1d<std::size_t> indices( nNode1D );
  for(std::size_t k=0; k<nNode1D; k++)
  {
    indices[k] = k * common.mStride + common.mIStart;
  }

  return indices;
}
