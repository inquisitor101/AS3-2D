#include "import_structure.hpp"




//-----------------------------------------------------------------------------------
// NImportFile namespace functions.
//-----------------------------------------------------------------------------------


void NImportFile::ImportAS3Grid
(
 CConfig   *config_container,
 CGeometry *geometry_container
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
			ImportAS3GridBinary(config_container, geometry_container); 
			break;
		}

		default: ERROR("Incorrect grid file format specified.");
	}

}

//-----------------------------------------------------------------------------------

void NImportFile::ImportAS3GridBinary
(
 CConfig   *config_container,
 CGeometry *geometry_container
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
	for( auto& zone: geometry_container->GetZoneGeometry() )
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
		EFaceElement imin = MapFaceElement.at( iconv[0] );
		EFaceElement imax = MapFaceElement.at( iconv[1] );
		EFaceElement jmin = MapFaceElement.at( iconv[2] );
		EFaceElement jmax = MapFaceElement.at( iconv[3] );

		// Ensure the correctness of the convention.
		if( (imin != EFaceElement::IMIN) || (imax != EFaceElement::IMAX) ||
				(jmin != EFaceElement::JMIN) || (jmax != EFaceElement::JMAX) )
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
		as3vector2d<unsigned int> buf_mark(nmark);
		as3vector2d<EFaceElement> buf_face(nmark);

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
		zone->InitializeMarkers(buf_mark, buf_face, buf_name);

		// Close the file.
		std::fclose(fh);
	}

	// Report output.
	std::cout << "Done." << std::endl;
}



