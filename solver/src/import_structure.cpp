#include "import_structure.hpp"

#include <set>


//-----------------------------------------------------------------------------------
// NImportFile namespace functions.
//-----------------------------------------------------------------------------------


void NImportFile::ImportAS3Grid
(
 CConfig            *config_container,
 CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that imports an AS3 grid file.
	*/
{
	// First, ensure this is indeed an AS3 file.
	if( config_container->GetMeshFormat() != EMeshFormat::AS3 )
	{
		ERROR("Grid must use an AS3 format.");
	}

	// Check which file format to use.
	switch( config_container->GetInputGridFormat() )
	{
		case(EFormatFile::BINARY):
		{
			ImportAS3GridBinary(config_container, multizone_geometry_container); 
			break;
		}

		default: ERROR("Incorrect grid file format specified.");
	}

}

//-----------------------------------------------------------------------------------

void NImportFile::ImportAS3GridBinary
(
 CConfig            *config_container,
 CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that imports an AS3 grid file in binary format.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Importing AS3 binary grid file in: " << std::endl;


	// For readability, abbreviate some variables.
	const int nstr = CGNS_STRING_SIZE;
	
	// Loop over the zones and import each grid separately.
	for( auto& zone: multizone_geometry_container->GetSinglezoneGeometry() )
	{
		// Get the filename in this zone.
		std::string filename = zone->GetGridFile();

		// Report output.
		std::cout << "  filename: " << filename << std::endl;

		// Standard error message.
		std::string errormsg = filename + " is not an AS3 binary file.";

		// Open the file for binary reading.
		FILE *fh = std::fopen( filename.c_str(), "rb");

		// Check file can be opened.
		if( !fh ) ERROR("File could not be opened for import.");

		/*
		 * Read the header information.
		 */

		// Read the magic number.
		unsigned int buf_uint;
		if( std::fread( &buf_uint, sizeof(unsigned int), 1, fh ) != 1 ) ERROR(errormsg);

		// Check if byte-swapping is needed.
		const bool swap = CheckByteSwapping(AS3_MAGIC_NUMBER, buf_uint);

		// Read the dimension of the problem, must be two.
		if( std::fread( &buf_uint, sizeof(unsigned int), 1, fh ) != 1 ) ERROR(errormsg);
		// Swap bytes, if necessary.
		if( swap ) NInputUtility::SwapBytes(&buf_uint, sizeof(unsigned int), 1);
		// Ensure this is a 2D grid file.
		if( buf_uint != 2 ) ERROR(filename + " must be 2D.");

		// Read the number of elements in the x and y-dimensions.
		unsigned int nelem[2];
		if( std::fread( &nelem, sizeof(unsigned int), 2, fh ) != 2 ) ERROR(errormsg);
		// Swap bytes, if necessary.
		if( swap ) NInputUtility::SwapBytes(&nelem, sizeof(unsigned int), 2);

		// Read the number of nodes2D in an element.
		unsigned int nnode2d;
		if( std::fread( &nnode2d, sizeof(unsigned int), 1, fh ) != 1 ) ERROR(errormsg);
		// Swap bytes, if necessary.
		if( swap ) NInputUtility::SwapBytes(&nnode2d, sizeof(unsigned int), 1);

		/*
		 * Read the coordinates.
		 */

		// Read the x-coordinates. Precision assumed double.
		size_t ntot = nelem[0]*nelem[1];
		as3vector2d<double> buf_x( ntot, as3vector1d<double>(nnode2d) );

		// Read x-coordinates in each element separately.
		for( auto& xx: buf_x )
		{
			if( std::fread( xx.data(), sizeof(double), nnode2d, fh ) != nnode2d ) ERROR(errormsg);
			// Swap bytes, if necessary.
			if( swap ) NInputUtility::SwapBytes(xx.data(), sizeof(double), nnode2d);
		}
	
		// Read the y-coordinates. Precision assumed double.
		as3vector2d<double> buf_y( ntot, as3vector1d<double>(nnode2d) );

		// Read y-coordinates in each element separately.
		for( auto& yy: buf_y )
		{
			if( std::fread( yy.data(), sizeof(double), nnode2d, fh ) != nnode2d ) ERROR(errormsg);
			// Swap bytes, if necessary.
			if( swap ) NInputUtility::SwapBytes(yy.data(), sizeof(double), nnode2d);
		}

		/*
		 * Read the markers.
		 */

		// Read the local element face index convention.
		unsigned int iconv[4];
		if( std::fread( &iconv, sizeof(unsigned int), 4, fh ) != 4 ) ERROR(errormsg);
		// Swap bytes, if necessary.
		if( swap ) NInputUtility::SwapBytes(&iconv, sizeof(unsigned int), 4);

		// Check that the map contains the indicial values.
		for(size_t i=0; i<4; i++)
		{
			if( !MapFaceElement.contains(iconv[i]) ) 
				ERROR("Face convention is wrong.");
		}

		// Explicitly map the values, according to the expected written convention.
		EFaceLocation imin = MapFaceElement.at( iconv[0] );
		EFaceLocation imax = MapFaceElement.at( iconv[1] );
		EFaceLocation jmin = MapFaceElement.at( iconv[2] );
		EFaceLocation jmax = MapFaceElement.at( iconv[3] );

		// Ensure the correctness of the convention.
		if( (imin != EFaceLocation::IMIN) || (imax != EFaceLocation::IMAX) ||
				(jmin != EFaceLocation::JMIN) || (jmax != EFaceLocation::JMAX) )
		{
			ERROR(filename + " adopts a difference element face convention.");
		}

		// Read the number of markers.
		unsigned int nmark;
		if( std::fread( &nmark, sizeof(unsigned int), 1, fh ) != 1 ) ERROR(errormsg);
		// Swap bytes, if necessary.
		if( swap ) NInputUtility::SwapBytes(&nmark, sizeof(unsigned int), 1);

		// Initialize vector of marker names.
		as3vector1d<std::string> buf_name;

		// Allocate vectors of marker indices and local face orientation.
		as3vector2d<unsigned int>  buf_mark(nmark);
		as3vector2d<EFaceLocation> buf_face(nmark);

		// Loop over each marker and read its information.
		for(size_t i=0; i<nmark; i++)
		{
			// Read the marker tag name. 
			char buff[nstr];
			if( std::fread( &buff[0], sizeof(char), nstr, fh ) != nstr ) ERROR(errormsg);
			//std::string name(buff);
			buf_name.push_back( buff );

			// Read the number of faces on this marker.
			unsigned int nf;
			if( std::fread( &nf, sizeof(unsigned int), 1, fh ) != 1 ) ERROR(errormsg);
			// Swap bytes, if necessary.
			if( swap ) NInputUtility::SwapBytes(&nf, sizeof(unsigned int), 1);

			// Resize the marker and copy the values to it.
			as3vector2d<unsigned int> imark( nf, as3vector1d<unsigned int>(2) );
			for( auto& m: imark )
			{
				if( std::fread( m.data(), sizeof(unsigned int), 2, fh ) != 2 ) ERROR(errormsg);
				// Swap bytes, if necessary.
				if( swap ) NInputUtility::SwapBytes(m.data(), sizeof(unsigned int), 2);
			}

			// Copy the values to their corresponding buffer.
			buf_mark[i].resize(nf);
			buf_face[i].resize(nf);
			for(size_t j=0; j<nf; j++)
			{
				buf_mark[i][j] = imark[j][0];
				if( !MapFaceElement.contains(imark[j][1]) ) ERROR("Face index could not be mapped.");
				buf_face[i][j] = MapFaceElement.at( imark[j][1] );
			}

		}

		// Initialize elements in this zone.
		zone->InitializeElements(buf_x, buf_y, nelem[0], nelem[1]);

		// Initialize markers in this zone.
		zone->InitializeMarkers(config_container, buf_mark, buf_face, buf_name);

		// Close the file.
		std::fclose(fh);
	}

	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

NImportFile::NAS3BinaryFile::CAS3BinaryMetadata NImportFile::NAS3BinaryFile::ReadAS3BinaryMetadata
(
 const std::string &filename
)
 /*
	* Function that imports an AS3 meta data in binary format.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Temporary lambdo to read one mapping, checking for duplicate IDs and names.
  auto lReadMapping = [&]() -> std::map<as3fileuint, std::string>
  {
    const auto count = ReadAS3UInt(file);

    std::map<as3fileuint, std::string> entries;
    std::set<std::string> names;

    for(as3fileuint i=0; i<count; i++)
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
  for(size_t i=0; i<nZone; i++)
  {
    // Read the relevant information.
    const size_t iZone  = ReadAS3UInt(file);
    const size_t nPoly  = ReadAS3UInt(file);
    const size_t niElem = ReadAS3UInt(file);
    const size_t njElem = ReadAS3UInt(file);

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
	* Function that imports an AS3 single zone grid in binary format.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Extract the header information.
  const size_t nPoly  = ReadAS3UInt(file);
  const size_t niElem = ReadAS3UInt(file);
  const size_t njElem = ReadAS3UInt(file);
  
  const bool isAffine = ReadBoolean(file);
 
  const size_t input_nodal_distribution = ReadAS3UInt(file);

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
	* Function that imports an AS3 external marker file in binary format.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Read the number of different marker families.
  const size_t nFamily = ReadAS3UInt(file);
 
  // Reserve the needed memory for the external marker family object.
  as3vector1d<CExternalFamilyMarker> external_family_marker;
  external_family_marker.reserve( nFamily ); 

  // Extract the relevant marker information for each family.
  for(size_t i=0; i<nFamily; i++)
  {
    const size_t iFamily   = ReadAS3UInt(file);
    const size_t nBoundary = ReadAS3UInt(file);
    const auto   name      = ReadString(file);

    // Consistency checks.
    if( iFamily != i ) ERROR(filename + " uses different family indices than expected.");

    // Initialize a family, and obtain a reference for it.
    auto& family = external_family_marker.emplace_back( std::move(name), nBoundary );

    // Get a reference for the markers, to initialize them later.
    auto& markers = family.GetMarkers();

    // Loop over each marker in this family.
    for(size_t j=0; j<nBoundary; j++)
    {
      const size_t iBoundary           = ReadAS3UInt(file);
      const size_t iZone               = ReadAS3UInt(file);
      const size_t input_face_location = ReadAS3UInt(file); 

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
	* Function that imports an AS3 internal marker file in binary format.
	*/
{
  // Open the AS3 file for binary reading.
  const auto file = OpenAS3BinaryFile(filename);

  // Read the number of internal markers, which are interfaces.
  const size_t nInterface = ReadAS3UInt(file);
 
  // This is always a single family, but we keep it as a vector, 
  // for general implementations.
  const size_t nFamily = 1;

  // Reserve the needed memory for the internal marker family object.
  as3vector1d<CInternalFamilyMarker> internal_family_marker;
  internal_family_marker.reserve( nFamily ); 

  // Extract the relevant marker information for each family.
  for(size_t i=0; i<nFamily; i++)
  {
    // Initialize a family, and obtain a reference for it.
    auto& family = internal_family_marker.emplace_back( nInterface );

    // Get a reference for the markers, to initialize them later.
    auto& markers = family.GetMarkers();

    // Loop over each marker in this family.
    for(size_t j=0; j<nInterface; j++)
    {
      const size_t iInterface           = ReadAS3UInt(file);
      
      const size_t iZone                = ReadAS3UInt(file);
      const size_t input_iface_location = ReadAS3UInt(file); 

      const size_t jZone                = ReadAS3UInt(file);
      const size_t input_jface_location = ReadAS3UInt(file); 

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
  *
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

NImportFile::NAS3BinaryFile::as3fileuint NImportFile::NAS3BinaryFile::ReadAS3UInt
(
 const CBinaryFile &file
)
 /*
  *
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
  return value;  
}

//-----------------------------------------------------------------------------------

std::string NImportFile::NAS3BinaryFile::ReadString
(
 const CBinaryFile &file
)
 /*
  *
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
  *
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
  * Infer coordinate precision from the remaining coordinate payload.
  * Restore the file position before returning.
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
