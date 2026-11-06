#include "geometry/CMultizoneGeometry.hpp"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "CConfig.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "io/NImport.hpp"
#include "NLogger.hpp"
#include "parallel/COpenMP.hpp"


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

  // Consistency check.
  if( mNZone != config_container->GetnZone() )
  {
    ERROR("Number of zones in imported grid does not match input by user in configuration.");
  }

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
