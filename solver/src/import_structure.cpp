#include "import_structure.hpp"

#include <cstddef>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include <set>


//-----------------------------------------------------------------------------------
// NImportFile namespace functions.
//-----------------------------------------------------------------------------------

NImportFile::NAS3BinaryFile::CAS3BinaryMetadata NImportFile::NAS3BinaryFile::ReadAS3BinaryMetadata
(
 const std::string &filename
)
 /*
	* Function that reads and validates the AS3 grid header and identifier mappings.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Temporary lambdo to read one mapping, checking for duplicate IDs and names.
  auto lReadMapping = [&]() -> std::map<std::size_t, std::string>
  {
    const auto count = ReadAS3UInt(file);

    std::map<std::size_t, std::string> entries;
    std::set<std::string> names;

    for(std::size_t i=0; i<count; i++)
    {
      const auto name = ReadString(file);
      const auto id   = ReadAS3UInt(file);
      
      if(!entries.emplace(id, name).second)
      {
        ERROR(filename + ": duplicate mapping ID.");
      } 
      
      if(!names.emplace(name).second)
      {
        ERROR(filename + ": duplicate mapping name.");
      }
    }

    return entries;
  }; 

  // Instantiate a metadata object.
  CAS3BinaryMetadata metadata;

  // Initialize the AS3 grid properties.
  metadata.mFilename = file.mFilename;
  metadata.mByteSwap = file.mSwap;
  metadata.mNDim     = ReadAS3UInt(file);

  // Consistency check.
  if(metadata.mNDim != 2) ERROR(filename + " must be 2D.");

  const std::size_t nZone = ReadAS3UInt(file);
  if(nZone == 0) ERROR(filename + " must contain at least one block.");

  // Reserve the number of single-zone grid metadata.
  metadata.mZones.reserve(nZone);

  // Read all per-zone metadata.
  for(std::size_t i=0; i<nZone; i++)
  {
    // Read the relevant information.
    const std::size_t iZone  = ReadAS3UInt(file);
    const std::size_t nPoly  = ReadAS3UInt(file);
    const std::size_t niElem = ReadAS3UInt(file);
    const std::size_t njElem = ReadAS3UInt(file);

    // Consistency checks.
    if( iZone != i ) ERROR(filename + " uses different zone indices than expected.");
    if( nPoly == 0 || niElem == 0 || njElem == 0 ) ERROR(filename + " has invalid block metadata."); 

    // Create a zone object with these information.
    CAS3ZoneMetadata zone{ iZone, nPoly, niElem, njElem };

    // Book-keep the information.
    metadata.mZones.push_back( zone ); 
  }


  // Extract a reference to the mappings object in the metadata class.
  auto& mappings = metadata.mMappings;

  // Read the edge mapping.
  for (const auto& entry : lReadMapping())
  {
    const auto id    = entry.first;
    const auto& name = entry.second;

    if      (name == "imin") mappings.mFaces.emplace(id, EFaceLocation::IMIN);
    else if (name == "imax") mappings.mFaces.emplace(id, EFaceLocation::IMAX);
    else if (name == "jmin") mappings.mFaces.emplace(id, EFaceLocation::JMIN);
    else if (name == "jmax") mappings.mFaces.emplace(id, EFaceLocation::JMAX);
    else if (name == "kmin" || name == "kmax") continue;
    else ERROR(filename + ": unknown face name: " + name);
  }

  // Consistency check.
  if(mappings.mFaces.size() != 4) ERROR(filename + ": missing required 2D face mappings.");

  // Read the boundary mapping.
  for( const auto& entry : lReadMapping() )
  {
    const auto id    = entry.first;
    const auto& name = entry.second;

    if      (name == "internal") mappings.mBoundaries.emplace(id, ETypeZoneMarker::INTERNAL);
    else if (name == "external") mappings.mBoundaries.emplace(id, ETypeZoneMarker::EXTERNAL);
    else if (name == "periodic") continue;
    else ERROR(filename + ": unknown boundary name: " + name);
  }

  // Consistency check.
  if(mappings.mBoundaries.size() != 2) ERROR(filename + ": missing required boundary mappings.");

  // Read the nodal distribution mapping.
  for (const auto& entry : lReadMapping())
  {
    const auto id    = entry.first;
    const auto& name = entry.second;

    const auto it = MapTypeDOF.find(name);

    if(it == MapTypeDOF.end()) ERROR(filename + ": unknown nodal distribution: " + name);

    mappings.mNodalDistributions.emplace(id, it->second);
  }

  if (mappings.mNodalDistributions.size() != MapTypeDOF.size())
  {
    ERROR(filename + ": missing required nodal distribution mappings.");
  }

  // Return the metadata object via copy elision.
  return metadata;
}

//-----------------------------------------------------------------------------------

std::unique_ptr<CSinglezoneGeometry> NImportFile::NAS3BinaryFile::ReadAS3BinarySinglezoneGrid
(
 const std::string        &filename,
 const CAS3BinaryMetadata &metadata,
 std::size_t               iZone
)
 /*
	* Function that imports one zone and its element coordinates from an AS3 binary file.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Extract the header information.
  const std::size_t nPoly  = ReadAS3UInt(file);
  const std::size_t niElem = ReadAS3UInt(file);
  const std::size_t njElem = ReadAS3UInt(file);
  
  const bool isAffine = ReadBoolean(file);
 
  const std::size_t input_nodal_distribution = ReadAS3UInt(file);

  // Deduce the type of nodal points.
  ETypeDOF nodal_distribution = metadata.mMappings.mNodalDistributions.at( input_nodal_distribution );

  // Extract the current zone's meta data.
  const auto& zone_metadata = metadata.mZones.at(iZone);

  // Check the block header against its meta data.
  if( nPoly  != zone_metadata.mNPoly  ||
      niElem != zone_metadata.mNiElem ||
      njElem != zone_metadata.mNjElem )
  {
    ERROR(filename + ": block header disagrees with metadata.");
  }
  

  // Construct the zone and allocate its final coordinate storage.
  auto grid = std::make_unique<CSinglezoneGeometry>(filename,
                                                    iZone,
                                                    nPoly,
                                                    niElem,
                                                    njElem,
                                                    nodal_distribution,
                                                    isAffine);


  // Determine the number of expected coordinates.
  const std::size_t nCoor = grid->GetnDim() * grid->GetnElem() * grid->GetnNodeGrid2D();

  // Deduce the coordinate precision
  const auto coor_bytes = DetermineCoordinatePrecisionBytes(file, nCoor);

  // Select the appropriate function to initialize the elements.
  switch( coor_bytes )
  {
    case(4): { ReadElementCoordinates<float> (file, *grid); break; }
    case(8): { ReadElementCoordinates<double>(file, *grid); break; }
    default: ERROR("Unsupported coordinate precision: " + filename);
  }

  // Return the grid via copy-elision.
  return grid;
}

//-----------------------------------------------------------------------------------

as3vector1d<CExternalFamilyMarker> NImportFile::NAS3BinaryFile::ReadAS3ExternalMarkers
(
 const std::string        &filename,
 const CAS3BinaryMetadata &metadata
)
 /*
	* Function that reads named physical boundary marker families from the AS3 grid file.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Read the number of different marker families.
  const std::size_t nFamily = ReadAS3UInt(file);
 
  // Reserve the needed memory for the external marker family object.
  as3vector1d<CExternalFamilyMarker> external_family_marker;
  external_family_marker.reserve( nFamily ); 

  // Extract the relevant marker information for each family.
  for(std::size_t i=0; i<nFamily; i++)
  {
    const std::size_t iFamily   = ReadAS3UInt(file);
    const std::size_t nBoundary = ReadAS3UInt(file);
    const auto   name      = ReadString(file);

    // Consistency checks.
    if( iFamily != i ) ERROR(filename + " uses different family indices than expected.");

    // Initialize a family, and obtain a reference for it.
    auto& family = external_family_marker.emplace_back( std::move(name), nBoundary );

    // Get a reference for the markers, to initialize them later.
    auto& markers = family.GetMarkers();

    // Loop over each marker in this family.
    for(std::size_t j=0; j<nBoundary; j++)
    {
      const std::size_t iBoundary           = ReadAS3UInt(file);
      const std::size_t iZone               = ReadAS3UInt(file);
      const std::size_t input_face_location = ReadAS3UInt(file); 

      // Consistency checks.
      if( iBoundary != j ) ERROR(filename + " uses different boundary indices than expected.");

      // Get the mapped face location.
      const EFaceLocation face_location = metadata.mMappings.mFaces.at( input_face_location ); 
    
      // Initialize a new marker with this information.
      markers.push_back( {iZone, face_location} );
    }
  }

  // Return the markers via copy-elision.
  return external_family_marker;
}

//-----------------------------------------------------------------------------------

as3vector1d<CInternalFamilyMarker> NImportFile::NAS3BinaryFile::ReadAS3InternalMarkers
(
 const std::string        &filename,
 const CAS3BinaryMetadata &metadata
)
 /*
	* Function that reads paired-zone interface marker families from the AS3 grid file.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Read the number of internal markers, which are interfaces.
  const std::size_t nInterface = ReadAS3UInt(file);
 
  // This is always a single family, but we keep it as a vector, 
  // for general implementations.
  const std::size_t nFamily = 1;

  // Reserve the needed memory for the internal marker family object.
  as3vector1d<CInternalFamilyMarker> internal_family_marker;
  internal_family_marker.reserve( nFamily ); 

  // Extract the relevant marker information for each family.
  for(std::size_t i=0; i<nFamily; i++)
  {
    // Initialize a family, and obtain a reference for it.
    auto& family = internal_family_marker.emplace_back( nInterface );

    // Get a reference for the markers, to initialize them later.
    auto& markers = family.GetMarkers();

    // Loop over each marker in this family.
    for(std::size_t j=0; j<nInterface; j++)
    {
      const std::size_t iInterface           = ReadAS3UInt(file);
      
      const std::size_t iZone                = ReadAS3UInt(file);
      const std::size_t input_iface_location = ReadAS3UInt(file); 

      const std::size_t jZone                = ReadAS3UInt(file);
      const std::size_t input_jface_location = ReadAS3UInt(file); 

      const bool is_reversed            = ReadBoolean(file);

      // Consistency checks.
      if( iInterface != j ) ERROR(filename + " uses different interface indices than expected.");

      // Get the mapped face locations.
      const EFaceLocation iface_location = metadata.mMappings.mFaces.at( input_iface_location ); 
      const EFaceLocation jface_location = metadata.mMappings.mFaces.at( input_jface_location ); 

      // Initialize a new marker with this information.
      markers.push_back( {iZone, jZone, iface_location, jface_location, is_reversed} );
    }
  }

  // Return the markers via copy-elision.
  return internal_family_marker;
}

//-----------------------------------------------------------------------------------

NImportFile::NAS3BinaryFile::CBinaryFile NImportFile::NAS3BinaryFile::OpenAS3BinaryFile
(
 const std::string &filename
)
 /*
  * Function that opens an AS3 binary file and determines its byte order.
  */
{
  // Create a binary file object.
  CBinaryFile file;

  // Assign its filename.
  file.mFilename = filename;

  // Assign its file handler too.
  file.mHandle.reset( std::fopen( filename.c_str(), "rb" ) );

  // Check if the file can be openned.
  if( !file.mHandle ) ERROR("Cannot open " + filename);

  // Read the magic number.
	as3fileuint magic{};
	if( std::fread( &magic, sizeof(magic), 1, file.mHandle.get() ) != 1 )
  {
    ERROR(filename + " is not an AS3 binary file.");
  }
  
  // Check if byte-swapping is needed.
  file.mSwap = NImportFile::CheckByteSwapping(AS3_MAGIC_NUMBER, magic);

  // Return the instance of this object.
  return file;
}

//-----------------------------------------------------------------------------------

std::size_t NImportFile::NAS3BinaryFile::ReadAS3UInt
(
 const CBinaryFile &file
)
 /*
  * Function that reads a 32-bit AS3 unsigned integer and converts it to host byte order.
  */
{
  as3fileuint value{};
  
  if( std::fread( &value, sizeof(value), 1, file.mHandle.get() ) != 1 )
  {
    ERROR(file.mFilename + " is not an AS3 binary file.");
  }
  
  if( file.mSwap )
  {
    NInputUtility::SwapBytes(&value, sizeof(value), 1);
  }
  return static_cast<std::size_t>(value);  
}

//-----------------------------------------------------------------------------------

std::string NImportFile::NAS3BinaryFile::ReadString
(
 const CBinaryFile &file
)
 /*
  * Function that reads a length-prefixed string from the AS3 binary file.
  */
{
  const auto length = ReadAS3UInt(file);
  std::string value(length, '\0');
  
  if(length != 0 && std::fread(&value[0], 1, value.size(), file.mHandle.get()) != value.size())
  {
    ERROR(file.mFilename + " is not an AS3 binary file.");
  }
  
  return value;
}

//-----------------------------------------------------------------------------------

bool NImportFile::NAS3BinaryFile::ReadBoolean
(
 const CBinaryFile &file
)
 /*
  * Function that reads a boolean value from the AS3 binary file.
  */
{
  const auto value = ReadAS3UInt(file);
  if( value > 1 ) ERROR("Invalid Boolean in " + file.mFilename);
  return static_cast<bool>(value);
}

//-----------------------------------------------------------------------------------

std::size_t NImportFile::NAS3BinaryFile::DetermineCoordinatePrecisionBytes
(
  const CBinaryFile &file,
  std::size_t        nCoor
)
 /*
  * Function that infers coordinate precision without changing the file position.
  */
{
  const auto& filename = file.mFilename;
  std::FILE* fh = file.mHandle.get();

  if( !fh ) ERROR("Invalid file handle: " + filename);
  if( nCoor == 0 ) ERROR("Invalid coordinate count: " + filename);

  const long headerEnd = std::ftell(fh);
  if( headerEnd < 0 ) ERROR("Cannot determine position in " + filename);

  if( std::fseek(fh, 0, SEEK_END) != 0 )
  {
    ERROR("Cannot seek to the end of " + filename);
  }

  const long fileEnd = std::ftell(fh);

  if( std::fseek(fh, headerEnd, SEEK_SET) != 0 )
  {
    ERROR("Cannot restore position in " + filename);
  }

  if( fileEnd < headerEnd ) ERROR("Invalid file length: " + filename);

  const auto remainingBytes =
    static_cast<std::size_t>(fileEnd - headerEnd);

  if( remainingBytes == 0 || remainingBytes % nCoor != 0 )
  {
    ERROR("Invalid coordinate payload size: " + filename);
  }

  return remainingBytes / nCoor;
}
