#include "geometry_structure.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <sstream>
#include <string>
#include "import_structure.hpp"
#include "log_structure.hpp"

#include <set>



//-----------------------------------------------------------------------------------
// CMultizoneGeometry member functions.
//-----------------------------------------------------------------------------------


CMultizoneGeometry::CMultizoneGeometry
(
 const CConfig *config_container
)
	:
		mNZone( config_container->GetnZone() )
 /*
	* Constructor that checks and imports the multizone grid.
	*/
{
	// Check and report output for the existance of the specified grid files.
	CheckExistanceGridFiles(config_container);

  // Import the grid.
  ImportGrid(config_container);

  // TODO: use the external and internal marker information to initialize the faces.
}

//-----------------------------------------------------------------------------------

CMultizoneGeometry::~CMultizoneGeometry
(
 void
)
 /*
	* Destructor for CMultizoneGeometry.
	*/
{

}

//-----------------------------------------------------------------------------------

void CMultizoneGeometry::ImportGrid
(
 const CConfig *config_container
)
 /*
  * Function that imports AS3 metadata, zone coordinates and zone-edge markers.
  */
{
	if( config_container->GetMeshFormat() != EMeshFormat::AS3 )
	{
		ERROR("Grid must use an AS3 format.");
	}

  if( config_container->GetInputGridFormat() != EFormatFile::BINARY )
  {
    ERROR("Currently, only binary AS3 grids can be imported.");
  }


  // For convenience, bring the import namespace in this scope.
  using namespace NImportFile::NAS3BinaryFile;

  // Extract the grid directory.
  std::string dir = config_container->GetGridDirectory();

  // Ensure the directory has the proper backslash.
  if (!dir.empty() && dir.back() != '/') dir += '/';

  // Obtain a metadata object for this grid.
  const auto metadata = ReadAS3BinaryMetadata(dir + "meta_data.bin");

  // Display information.
  NLogger::DisplayAS3BinaryMetadata(metadata, std::cout);

  // Deduce the number of zones.
  mNZone = metadata.mZones.size();

  // Reserve the correct amount of memory for the grid.
  mSinglezoneGeometry.reserve( mNZone );

  // Obtain the single zone grids.
  for(std::size_t iZone=0; iZone<mNZone; iZone++)
  {
    const auto& filename = dir + "blocks/block_" + std::to_string(iZone) + ".bin"; 
    mSinglezoneGeometry.push_back( ReadAS3BinarySinglezoneGrid(filename, metadata, iZone) );
  }
 

  // Import the internal markers.
  mInternalFamilyMarkers = ReadAS3InternalMarkers(dir + "interface_boundaries.bin", metadata);

  // Import the external markers.
  mExternalFamilyMarkers = ReadAS3ExternalMarkers(dir + "external_boundaries.bin", metadata);


  // TODO: The below should be in a function here. 

  // Process the periodic boundaries, specified by the user and remove them from external markers. 
  const auto& periodic_param = config_container->GetPeriodicParamMarker(); 

  const std::size_t nPeriodic = periodic_param.size();

  mPeriodicFamilyMarkers.reserve( nPeriodic );

  // Search for an external family by name.
  auto FindExternalFamily = [&](const std::string &name)
  {
    return std::find_if(
      mExternalFamilyMarkers.begin(),
      mExternalFamilyMarkers.end(),
      [&name](const CExternalFamilyMarker &family)
      {
        return family.GetName() == name;
      }
    );
  };

  for( const auto& marker : periodic_param )
  {
    const auto& iname = marker->mName;
    const auto& jname = marker->mNameMatching;
  
    // Each periodic pair must refer to two different families.
    if( iname == jname )
    {
      ERROR("Periodic boundaries must have different family names: " + iname);
    }
  
    // Find the I family.
    const auto itI = FindExternalFamily(iname);
    if( itI == mExternalFamilyMarkers.end() )
    {
      ERROR("Cannot find periodic boundary family: " + iname);
    }
  
    // Find the matching J family.
    const auto itJ = FindExternalFamily(jname);
    if( itJ == mExternalFamilyMarkers.end() )
    {
      ERROR("Cannot find periodic boundary family: " + jname);
    }
  
    // Construct the periodic family before invalidating either iterator.
    mPeriodicFamilyMarkers.emplace_back(this, *itI, *itJ, marker->mVectorTrans);
  
    // Erase the later element first, preserving the earlier iterator.
    if( itI < itJ )
    {
      mExternalFamilyMarkers.erase(itJ);
      mExternalFamilyMarkers.erase(itI);
    }
    else
    {
      mExternalFamilyMarkers.erase(itI);
      mExternalFamilyMarkers.erase(itJ);
    }
  }


  // Display information.
  NLogger::DisplayInternalMarkers(mInternalFamilyMarkers);
  NLogger::DisplayExternalMarkers(mExternalFamilyMarkers);
  NLogger::DisplayPeriodicMarkers(mPeriodicFamilyMarkers); 
}

//-----------------------------------------------------------------------------------

void CMultizoneGeometry::CheckExistanceGridFiles
(
 const CConfig *config_container
)
 /*
	* Function that checks the configured grid directory and required AS3 files.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Generating geometry: " << std::endl;

	// Report output.
	std::cout << "  expecting " << mNZone << " grid files. " << std::endl;
	std::cout << "   ... detected: " << std::endl;

  std::string dir = config_container->GetGridDirectory();

	// Loop over all the expected zones and report their grid filenames.
	for(std::size_t i=0; i<mNZone; i++)
	{
		std::string filename = dir + "blocks/block_" + std::to_string(i) + ".bin";
		// Check if the file exists.
		std::ifstream file(filename);
		if( !file.good() )
		{
			std::ostringstream message;
			message << "Could not open file: "
				      << "'" << filename << "'";
			ERROR(message.str());
		}

		// Report output.
		std::cout << "     zone: " << i << "): " << filename << "\n";
	}

  // TODO: also report the internal and external boundaries and metadata.

	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void CMultizoneGeometry::InitializeGridTopology
(
 const COpenMP *openmp_container
)
 /*
	* Function that builds flattened element and face indices and load-balancing permutations.
	*/
{
	// Initialize the faces in the i-direction.
	mMultizoneIFaces.InitializeFaces(this);

	// Initialize the faces in the j-direction.
	mMultizoneJFaces.InitializeFaces(this);


  // Temporary lambda to flatten the face indices for a given multizone face container.
  const auto lFlattenFaceIndices = [](const auto& multizone_faces, auto& flattened_indices)
  {
    // Extract the number of faces.
    const auto nFaces = multizone_faces.GetnFacesTotal();
    
    // Reserve memory preciselty as needed.
    flattened_indices.reserve(nFaces);
 
    // Loop over each face and construct its indicial mapping.
    for(std::size_t i=0; i<nFaces; i++)
    {
      // Resolve the face kind, batch index and batch-local face index.
      flattened_indices.push_back( multizone_faces.ResolveFace(i) );
    }

  };


  // Initialize the flattened facial indices in the i- and j-directions.
  lFlattenFaceIndices(mMultizoneIFaces, mFlattenedIndexIFace);
  lFlattenFaceIndices(mMultizoneJFaces, mFlattenedIndexJFace);



	// Deduce the total number of elements in the entire multizone grid.
	mNElemTotal = 0;
	for( const auto& zone : mSinglezoneGeometry ) mNElemTotal += zone->GetnElem();

	// Reserve the needed memory for the element indices.
	mFlattenedIndexVolumeElement.reserve( mNElemTotal );

	// Deduce the flattening strategy of the elements.
	for(std::size_t iZone=0; iZone<mSinglezoneGeometry.size(); iZone++)
	{
		for(std::size_t iElem=0; iElem<mSinglezoneGeometry[iZone]->GetnElem(); iElem++)
		{
			mFlattenedIndexVolumeElement.emplace_back( CFlattenedElementIndex{iZone, iElem} ); 
		}
	}


  // For load-balancing reasons, select the number of chunks to be equal to the number of OpenMP threads. 
  const std::size_t nChunk = openmp_container->GetnThread();

  // Obtain a better load-balanced face partitioning strategy for the i-faces.
  mIFaceLoadBalancedPermutation = CLoadBalancedFacePermutation( mFlattenedIndexIFace, 
                                                                nChunk, 
                                                                EFaceLoadBalanceStrategy::GREEDY);
  
  // Obtain a better load-balanced face partitioning strategy for the j-faces.
  mJFaceLoadBalancedPermutation = CLoadBalancedFacePermutation( mFlattenedIndexJFace,
                                                                nChunk,
                                                                EFaceLoadBalanceStrategy::GREEDY);
}



//-----------------------------------------------------------------------------------
// CSinglezoneGeometry member functions.
//-----------------------------------------------------------------------------------


CSinglezoneGeometry::CSinglezoneGeometry
(
  const std::string &gridfile,
  std::size_t        iZone,
  std::size_t        nPoly,
  std::size_t        niElem,
  std::size_t        njElem,
  ETypeDOF           nodalDistribution,
  bool               isAffine
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

  const std::size_t nNode1D = std::size_t(mNPolyGrid) + 1;
  const std::size_t nNode2D = nNode1D * nNode1D;
  const std::size_t nElem = std::size_t(mNiElem) * mNjElem;

  GenerateNodalFaceIndices(); // TODO: make it explicitly depend on nPoly as input

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

void CSinglezoneGeometry::GenerateNodalFaceIndices
(
 void
)
 /*
	* Function that builds reference-node indices for the four local quadrilateral edges.
	*/
{
	// Deduce the number of solution points in 1D in this zone.
	const std::size_t nSol1D = mNPolyGrid+1;

	// TODO: change to a contiguous array, using CMatrixAS3.
  // Allocate memory for the nodal indices on all sides.
	mFaceNodalIndices.resize( 4, as3vector1d<std::size_t>(nSol1D) );

	// Offsets in the i- and j-direction.
	const std::size_t di = mNPolyGrid;
	const std::size_t dj = nSol1D*(nSol1D-1);
	
	// Loop over the face nodes and generate their local indices.
	for(std::size_t i=0; i<nSol1D; i++)
	{
		// Face: IMIN.
		mFaceNodalIndices[0][i] = i*nSol1D;
		// Face: IMAX.
		mFaceNodalIndices[1][i] = mFaceNodalIndices[0][i] + di;
		// Face: JMIN.
		mFaceNodalIndices[2][i] = i;
		// Face: JMAX.
		mFaceNodalIndices[3][i] = mFaceNodalIndices[2][i] + dj;
	}
}

//-----------------------------------------------------------------------------------

void CSinglezoneGeometry::InitializeElements
(
 as3vector2d<double> &xcoor,
 as3vector2d<double> &ycoor,
 std::size_t         niElem,
 std::size_t         njElem
)
 /*
	* Function that constructs element geometries from per-element coordinate vectors.
	*/
{
	// Deduce the number of elements and ensure consistency.
	if( xcoor.size() != ycoor.size() ) ERROR("Number of coordinates in x and y is not identical.");

	// Ensure the total number of elements is correct.
	if( niElem*njElem != xcoor.size() ) ERROR("Inconsistency in number of elements.");

	// Set the number of elements in each dimension.
	mNiElem = niElem;
	mNjElem = njElem;

	// Allocate memory for the total elements.
	mElementGeometry.resize( xcoor.size() );

	// Total number of grid nodes in 2D, as user-specified.
	const std::size_t nNode2D = (mNPolyGrid+1)*(mNPolyGrid+1);

	// Loop over each element and instantiate its coordinates.
	for(std::size_t i=0; i<xcoor.size(); i++)
	{
		// Ensure the polynomial order is correct.
		if( (xcoor[i].size() != nNode2D) || (ycoor[i].size() != nNode2D) )
    {
      ERROR("Elements do not match polynomial order.");
    }

		// Instantiate the current element.
		mElementGeometry[i] = std::make_unique<CElementGeometry>( xcoor[i], ycoor[i] );
	}
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


//-----------------------------------------------------------------------------------
// CElementGeometry member functions.
//-----------------------------------------------------------------------------------

CElementGeometry::CElementGeometry
(
 std::size_t nNode2D
)
  : mCoordSolDOFs(2, nNode2D)
 /*
  * Constructor that allocates geometry-node coordinates for one element.
  */
{

}

//-----------------------------------------------------------------------------------

CElementGeometry::CElementGeometry
(
 as3vector1d<double> &x,
 as3vector1d<double> &y
)
 /*
	* Constructor that copies the geometry-node coordinates.
	*/
{
	// Ensure consistency.
	if( x.size() != y.size() ) ERROR("Inconsistent size in coordinates.");

	// Allocate memory for the coordinates.
	mCoordSolDOFs.resize( 2, x.size() );

	// Copy the coordinates.
	for(std::size_t i=0; i<x.size(); i++)
	{
		mCoordSolDOFs(0,i) = x[i];
		mCoordSolDOFs(1,i) = y[i];
	}
}

//-----------------------------------------------------------------------------------

CElementGeometry::~CElementGeometry
(
 void
)
 /*
	* Destructor for CElementGeometry.
	*/
{

}


