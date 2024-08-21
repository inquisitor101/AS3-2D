#include "config_structure.hpp"


//-----------------------------------------------------------------------------------
// CConfig member functions.
//-----------------------------------------------------------------------------------


CConfig::CConfig
(
 const char *filename
)
 /*
	* Constructor for the configuration class, which defines the simulation parameters.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Reading configuration file: " << filename << std::endl;

	// Message stream.
	std::ostringstream message;

  // Check if file exists.
  std::ifstream inputfile(filename);
  if( !inputfile.good() ) ERROR("Configuration file could not be opened!");
	
	// Close the file, since it is being opened at each subsequent step.
  inputfile.close();


	// Read the zone connectivity options.
	if( !ReadZoneConnectivityOptions(filename) )
	{
		message << "Failed to extract zone connectivity options from " << filename;
		ERROR(message.str());
	}

	// Read the zone configuration options.
	if( !ReadZoneConfigurationOptions(filename) )
	{
		message << "Failed to extract zone configuration options from " << filename;
		ERROR(message.str());
	}

	// Read the output options.
	if( !ReadOutputOptions(filename) )
	{
		message << "Failed to extract output options from " << filename;
		ERROR(message.str());
	}
	
	// Read the temporal options.
	if( !ReadTemporalOptions(filename) )
	{
		message << "Failed to extract temporal options from " << filename;
		ERROR(message.str());
	}

	// Read the initial condition options.
	if( !ReadInitialConditionOptions(filename) )
	{
		message << "Failed to extract initial condition options from " << filename;
		ERROR(message.str());
	}
	
	// Read the boundary condition options.
	if( !ReadBoundaryConditionOptions(filename) )
	{
		message << "Failed to extract boundary condition options from " << filename;
		ERROR(message.str());
	}
	


	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

CConfig::~CConfig
(
 void
)
 /*
	* Destructor, which cleans up after the configuration class.
	*/
{

}

//-----------------------------------------------------------------------------------

bool CConfig::ReadZoneConnectivityOptions
(
 const char *filename
)
 /*
	* Function that reads the zone connectivity options.
	*/
{
	// Report output.
	std::cout << "  reading zone connectivity............... ";	

  // Open input file.
  std::ifstream paramfile(filename);

  // Buffer for storing temporary strings.
	std::string buffer;

	// Read the mesh format specified.
	NInputUtility::AddScalarOption(paramfile, "MESH_FORMAT", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mMeshFormat = GenericScalarMap(MapMeshFormat, buffer, "MESH_FORMAT");

	// Read the total number of zones expected.
	NInputUtility::AddScalarOption(paramfile, "NUMBER_ZONE", mNZone, true);

	// Read the format of the input grid files.
	NInputUtility::AddScalarOption(paramfile, "INPUT_GRID_FORMAT", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mInputGridFormat = GenericScalarMap(MapFormatFile, buffer, "INPUT_GRID_FORMAT");

	// Close the file, since it will be opened later.
	paramfile.close();

	// Extract the information in the zone connectivity file.
	ExtractZoneGridFiles(filename);

	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

bool CConfig::ReadZoneConfigurationOptions
(
 const char *filename
)
 /*
	* Function that reads the zone configuration options.
	*/
{
	// Report output.
	std::cout << "  reading zone configurations............. ";

  // Open input file.
  std::ifstream paramfile(filename);

  // Buffer for storing temporary strings.
	as3vector1d<std::string> buffer;

  // Read the zone markers names.
	NInputUtility::AddVectorOption(paramfile, "MARKER_ZONE", buffer, true);
  // Pad entries for the buffer vector. This serves as an error check, rather than padding.
  PadEntriesVectorData(buffer, "MARKER_ZONE", mNZone);
	// Deduce the zone indices.
	NInputUtility::ExtractLastDigitFromString(buffer, mZoneIndex);

	// Read the type of solver in each zone.
	NInputUtility::AddVectorOption(paramfile, "TYPE_SOLVER", buffer, true);
  // Pad entries for the buffer vector.
  PadEntriesVectorData(buffer, "TYPE_SOLVER", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeSolver = GenericVectorMap(MapTypeSolver, buffer, "TYPE_SOLVER");

	// Default buffer layer value is: none.
	as3vector1d<std::string> dvalue = { "NONE" }; 
	// Read the type of buffer zone in each zone. Default: none.
	NInputUtility::AddVectorOption(paramfile, "TYPE_BUFFER_LAYER", buffer, dvalue, true);
  // Pad entries for the buffer vector.
  PadEntriesVectorData(buffer, "TYPE_BUFFER_LAYER", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeBufferLayer = GenericVectorMap(MapTypeBufferLayer, buffer, "TYPE_BUFFER_LAYER");

	// Read the type of DOFs in each zone.
	NInputUtility::AddVectorOption(paramfile, "TYPE_DOF", buffer, true);
  // Pad entries for the buffer vector.
  PadEntriesVectorData(buffer, "TYPE_DOF", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeDOF = GenericVectorMap(MapTypeDOF, buffer, "TYPE_DOF");

	// Read the grid polynomial orders in each zone.
	NInputUtility::AddVectorOption(paramfile, "POLY_ORDER_GRID", mNPolyGrid, true);
  // Pad entries for the mNPolyGrid vector.
  PadEntriesVectorData(mNPolyGrid, "POLY_ORDER_GRID", mNZone, 1, 2);

	// Read the solution polynomial orders in each zone.
	NInputUtility::AddVectorOption(paramfile, "POLY_ORDER_SOL", mNPolySol, true);
  // Pad entries for the mNPolySol vector.
  PadEntriesVectorData(mNPolySol, "POLY_ORDER_SOL", mNZone, 1, 2);

	// Consistency checks.
	ConsistencyCheckZoneConfiguration();

	// Close file.
  paramfile.close();

	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

bool CConfig::ReadOutputOptions
(
 const char *filename
)
 /*
	* Function that reads the output options.
	*/
{
	// Report output.
	std::cout << "  reading output information.............. ";

  // Open input file.
  std::ifstream paramfile(filename);

  // Buffers for storing temporary strings.
	std::string              buffer;
	as3vector1d<std::string> buffervec;


	// Read the output visualization format specified.
	NInputUtility::AddScalarOption(paramfile, "OUTPUT_VIS_FORMAT", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mOutputVisFormat = GenericScalarMap(MapVisualFormat, buffer, "OUTPUT_VIS_FORMAT");

	// Read the output visualization filename.
	NInputUtility::AddScalarOption(paramfile, "OUTPUT_VIS_FILENAME", mOutputVisFilename, true);
	// Read the output solution filename.
	NInputUtility::AddScalarOption(paramfile, "OUTPUT_SOL_FILENAME", mOutputSolFilename, true);

	// Reada the output visualization variables written.
	NInputUtility::AddVectorOption(paramfile, "WRITE_VIS_VARIABLE", buffervec, true);
	// Deduce the type from the buffer, based on the mapping.
	mWriteVisVar = GenericVectorMap(MapWriteVariable, buffervec, "WRITE_VIS_VARIABLE");

	// Close file.
  paramfile.close();

	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

bool CConfig::ReadTemporalOptions
(
 const char *filename
)
 /*
	* Function that reads the temporal options.
	*/
{
	// Report output.
	std::cout << "  reading temporal information............ ";
  
	// Open input file.
  std::ifstream paramfile(filename);

  // Buffer for storing temporary strings.
	std::string buffer;

	// Read the temporal discretization specified.
	NInputUtility::AddScalarOption(paramfile, "TEMPORAL_SCHEME", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mTemporalScheme = GenericScalarMap(MapTemporalScheme, buffer, "TEMPORAL_SCHEME");

	// Read the simulation starting time.
	NInputUtility::AddScalarOption(paramfile, "START_TIME", mStartTime, true);
	// Read the simulation ending time.
	NInputUtility::AddScalarOption(paramfile, "FINAL_TIME", mFinalTime, true);
	// Read the temporal step specified.
	NInputUtility::AddScalarOption(paramfile, "TIME_STEP", mTimeStep, true);
	// Read the maximum number of temporal iterations.
	NInputUtility::AddScalarOption(paramfile, "MAX_ITER_TIME", mMaxIterTime, true);

	// Close file.
  paramfile.close();

	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

bool CConfig::ReadInitialConditionOptions
(
 const char *filename
)
 /*
	* Function that reads the initial condition options.
	*/
{
	// Report output.
	std::cout << "  reading initial condition information... ";
  
	// Open input file.
  std::ifstream paramfile(filename);

  // Buffer for storing temporary strings.
	std::string buffer;

	// Read the type of initial condition specified.
	NInputUtility::AddScalarOption(paramfile, "TYPE_IC", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mTypeIC = GenericScalarMap(MapTypeIC, buffer, "TYPE_IC");

	// Close the file, since it is being opened during the processing of the ICs.
  paramfile.close();

	// Check what else to extract, based on the type of IC.
	switch(mTypeIC)
	{
		case(ETypeIC::GAUSSIAN_PRESSURE): { IC_GaussianPressure(filename); break; }
		case(ETypeIC::ISENTROPIC_VORTEX): { IC_IsentropicVortex(filename); break; }
		default: ERROR("Cannot extract additional information from the specified IC.");
	}


	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

bool CConfig::ReadBoundaryConditionOptions
(
 const char *filename
)
 /*
	* Function that reads the boundary condition options.
	*/
{
	// Report output.
  std::cout << "  reading boundary condition information.. ";
	
	// Check for and extract periodic boundary condition information, if present.
	ExtractInfoPeriodicBC(filename);


	// Check that all the markers are indeed unique. 
	if( mBoundaryMarkers.size() )
	{
		// Loop over the elements and ensure uniqueness.
		for(size_t i=0; i<mBoundaryMarkers.size(); i++)
			for(size_t j=i+1; j<mBoundaryMarkers.size(); j++)
				if( mBoundaryMarkers[i].first == mBoundaryMarkers[j].first ) ERROR("Markers are not unique."); 
	}
	else
	{
		ERROR("No boundary conditions have been deteced.");
	}


	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

void CConfig::ExtractInfoPeriodicBC
(
 const char *filename
)
 /*
	* Function that reads the periodic BC information, if present.
	*/
{
	// Open input file.
  std::ifstream paramfile(filename);

	// Default value, needed to check if/when this BC is not found in the file.
	as3vector1d<std::string> defval = { "ERROR404" };

	// Condition for when to stop searching for periodic BCs.
	bool condition = true;

	// Loop over the file and search for any periodic BC markers.
	while(condition)
	{
		// Buffers for storing temporary strings.
		as3vector1d<std::string> buffer;
		
		// Read the information of the periodic marker, if specified.
		NInputUtility::AddVectorOption(paramfile, "MARKER_BC_PERIODIC", buffer, defval, false);

		// Check if any new values are found.
		if( buffer != defval )
		{
			// Ensure the values are a pair.
			if( buffer.size() != 2 ) ERROR("Periodic BCs must specify two markers.");

			// Accumulate the information: iMarker, jMarker.
			mMarkerNamePeriodic.push_back( buffer );

			// The periodic markers are a pair: iMarker, jMarker. Also, include their reverse 
			// order: jMarker, iMarker. This is not strictly needed, but helps in the processing.
			std::reverse( buffer.begin(), buffer.end() );
			// Include the reversed order: jMarker, iMarker.
			mMarkerNamePeriodic.push_back( buffer );
		}
		else
		{
			// Prepare to exit from this loop.
			condition = false;
		}
	}


	// Check that all the markers are indeed unique. This can be done more elegantly, but this
	// is fine for now, at the preprocessing stage.
	if( mMarkerNamePeriodic.size() )
	{
		// Abbreviation, for readability.
		auto& m = mMarkerNamePeriodic;
		
		// Loop over the elements and ensure uniqueness.
		for(size_t i=0; i<m.size(); i++)
			for(size_t j=i+1; j<m.size(); j++)
				if( (m[i][0] == m[j][0]) || (m[i][1] == m[j][1]) ) ERROR("Periodic markers are not unique."); 
	}

	// Accumulate the marker name and type for book-keeping purposes.
	for(size_t i=0; i<mMarkerNamePeriodic.size(); i++)
	{
		mBoundaryMarkers.push_back( {mMarkerNamePeriodic[i][0], ETypeBC::PERIODIC} );
	}


	// Close file.
  paramfile.close();
}

//-----------------------------------------------------------------------------------

void CConfig::IC_GaussianPressure
(
 const char *filename
)
 /*
	* Function that extract the additional parameters in a Gaussian pressure pulse IC.
	*/
{
	// Open input file.
  std::ifstream paramfile(filename);

	// Values that store the respective information.
	as3double vratio;
	as3double vwidth;	
	as3double vmach;
	as3double vangle;
	as3vector1d<as3double> vcenter;

	// Read the pulse ratio specified.
	NInputUtility::AddScalarOption(paramfile, "DISTURBANCE_RATIO", vratio, true);
	// Read the pulse width specified.
	NInputUtility::AddScalarOption(paramfile, "DISTURBANCE_WIDTH", vwidth, true);
	
	// Read the coordinates of the center of the pulse.
	NInputUtility::AddVectorOption(paramfile, "DISTURBANCE_CENTER", vcenter, true);

	// Read the freestream Mach number specified.
	NInputUtility::AddScalarOption(paramfile, "FREESTREAM_MACH", vmach, true);
	// Read the freestream flow angle specified.
	NInputUtility::AddScalarOption(paramfile, "FREESTREAM_FLOW_ANGLE", vangle, true);

	// Ensure the center provided is two-dimensional.
	if( vcenter.size() != 2 ) ERROR("Pulse center must have x and y only.");

	// Assemble the vector containing the IC data.
	mDataIC.resize(6);

	// Pack the data into the vector. Use the same convention to unpack.
	mDataIC[0] = vratio;
	mDataIC[1] = vwidth;
	mDataIC[2] = vmach;
	mDataIC[3] = vangle;
	mDataIC[4] = vcenter[0];
	mDataIC[5] = vcenter[1];

	// Close file.
  paramfile.close();
}

//-----------------------------------------------------------------------------------

void CConfig::IC_IsentropicVortex
(
 const char *filename
)
 /*
	* Function that extract the additional parameters in an isentropic vortex IC.
	*/
{
	// Open input file.
  std::ifstream paramfile(filename);

	// Values that store the respective information.
	as3double vratio;
	as3double vwidth;	
	as3double vmach;
	as3double vangle;
	as3vector1d<as3double> vcenter;

	// Read the pulse ratio specified.
	NInputUtility::AddScalarOption(paramfile, "DISTURBANCE_RATIO", vratio, true);
	// Read the pulse width specified.
	NInputUtility::AddScalarOption(paramfile, "DISTURBANCE_WIDTH", vwidth, true);
	
	// Read the coordinates of the center of the pulse.
	NInputUtility::AddVectorOption(paramfile, "DISTURBANCE_CENTER", vcenter, true);

	// Read the freestream Mach number specified.
	NInputUtility::AddScalarOption(paramfile, "FREESTREAM_MACH", vmach, true);
	// Read the freestream flow angle specified.
	NInputUtility::AddScalarOption(paramfile, "FREESTREAM_FLOW_ANGLE", vangle, true);

	// Ensure the center provided is two-dimensional.
	if( vcenter.size() != 2 ) ERROR("Pulse center must have x and y only.");

	// Assemble the vector containing the IC data.
	mDataIC.resize(6);

	// Pack the data into the vector. Use the same convention to unpack.
	mDataIC[0] = vratio;
	mDataIC[1] = vwidth;
	mDataIC[2] = vmach;
	mDataIC[3] = vangle;
	mDataIC[4] = vcenter[0];
	mDataIC[5] = vcenter[1];

	// Close file.
  paramfile.close();
}

//-----------------------------------------------------------------------------------

void CConfig::ConsistencyCheckZoneConfiguration
(
 void
)
 /*
	* Function that does a consistency check for zone configuration.
	*/
{
	// Check that the zone indices are sequential.
	as3vector1d<unsigned short> tmp = mZoneIndex;
	// Sort the data sequentially.
	std::sort( tmp.begin(), tmp.end() );

	// Check that the data start from 0 and end at nZone-1.
	unsigned short j=0;
	for(unsigned short i=0; i<mNZone; i++)
	{
		if( tmp[i] != j++ ) ERROR("Zone indices must be sequential, starting from 0.");
	}				
}

//-----------------------------------------------------------------------------------

void CConfig::ExtractZoneGridFiles
(
 const char *filename
)
 /*
	* Function that extracts the grid zone connectivity information.
	*/
{
	// Message stream.
	std::ostringstream message;

  // Check if file exists.
  std::ifstream paramfile(filename);

	// Reserve memory for the grid filenames.
	mZoneGridFilename.resize(mNZone);

	// Loop over all the expected zones and read their grid filenames.
	for(unsigned short i=0; i<mNZone; i++)
	{
		// Current zone filename.
		std::string ifile = "GRID_FILENAME_ZONE_" + std::to_string(i);
		NInputUtility::AddScalarOption(paramfile, ifile.c_str(), mZoneGridFilename[i], true);

		// Ensure that the grid extension is correct.
		switch( mMeshFormat )
		{
			// Ensure a Plot3D grid is supplied.
			case( EMeshFormat::PLOT3D ): 
			{
				if( NInputUtility::GetFileExtension(mZoneGridFilename[i]) != "xyz" ) 
				{
					ERROR("Wrong mesh format detected. AS3 expects a Plot3D grid."); 
				}

				// This is not supported for now, issue an error.
				ERROR("PLOT3D is not (yet) supported.");
				break;
			}

			// Ensure an AS3 grid is supplied.
			case( EMeshFormat::AS3 ): 
			{
				if( NInputUtility::GetFileExtension(mZoneGridFilename[i]) != "as3" ) 
				{
					ERROR("Wrong mesh format detected. AS3 expects a native (AS3) grid."); 
				}
				break;
			}

			// Issue an error if format is unknown.
			default: ERROR("Unknown mesh format found.");
		}

		// Check if the file exists.
		std::ifstream file(mZoneGridFilename[i]);
		if( !file.good() )
		{
			std::ostringstream message;
			message << "Could not open file: "
				      << "'" << mZoneGridFilename[i] << "'";
			ERROR(message.str());
		}
	}

	// Close file.
  paramfile.close();
}

//-----------------------------------------------------------------------------------

ETypeBC CConfig::DetermineMarkerBC
(
 std::string &marker
)
 /*
	* Function that determines the boundary condition for a given marker.
	*/
{
	// Search for the relevant marker name.
	for( auto& bm: mBoundaryMarkers )
		if( bm.first == marker ) return bm.second;

	// If the code reaches this far, the marker has not been found.
	ERROR("Could not detect marker: " + marker);

	// The program should not reach this far. Return something to avoid a compiler error.
	return ETypeBC::PERIODIC;
}

//-----------------------------------------------------------------------------------
	
std::string CConfig::FindMatchingPeriodicMarker
(
 const std::string &imarker
)
 /*
	* Function that finds the name of the matching periodic boundary marker.
	*/
{
	// Search for the relevant marker name.
	for( auto& mp: mMarkerNamePeriodic )
		if( mp[0] == imarker ) return mp[1];

	// If the code reaches this far, the marker has not been found.
	ERROR("Could not find the matching marker for: " + imarker);

	// The program should not reach this far. Return something to avoid a compiler error.
	return "";
}





