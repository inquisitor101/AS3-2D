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
	std::string              buffer;
	as3vector1d<std::string> buffervec;

	// Read the total number of zones expected.
	NInputUtility::AddScalarOption(paramfile, "NUMBER_ZONE", mNZone, true);
 
	// Read the type of solver in each zone.
	NInputUtility::AddVectorOption(paramfile, "TYPE_SOLVER", buffervec, true);
  // Pad entries for the buffer vector.
  PadEntriesVectorData(buffervec, "TYPE_SOLVER", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeSolver = GenericVectorMap(MapTypeSolver, buffervec, "TYPE_SOLVER");

	// Default buffer layer value is: none.
	as3vector1d<std::string> dvalue = { "NONE" }; 
	// Read the type of buffer zone in each zone. Default: none.
	NInputUtility::AddVectorOption(paramfile, "TYPE_BUFFER_LAYER", buffervec, dvalue, true);
  // Pad entries for the buffer vector.
  PadEntriesVectorData(buffervec, "TYPE_BUFFER_LAYER", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeBufferLayer = GenericVectorMap(MapTypeBufferLayer, buffervec, "TYPE_BUFFER_LAYER");

	// Read the type of DOFs in each zone.
	NInputUtility::AddVectorOption(paramfile, "TYPE_DOF", buffervec, true);
  // Pad entries for the buffer vector.
  PadEntriesVectorData(buffervec, "TYPE_DOF", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeDOF = GenericVectorMap(MapTypeDOF, buffervec, "TYPE_DOF");

	// Read the solution/grid polynomial orders in each zone.
	NInputUtility::AddVectorOption(paramfile, "POLY_ORDER", mNPoly, true);
  // Pad entries for the mNPoly vector.
  PadEntriesVectorData(mNPoly, "POLY_ORDER", mNZone, 1, 2);

	// Read the type of Riemann solver specified in each zone.
	NInputUtility::AddVectorOption(paramfile, "RIEMANN_SOLVER", buffervec, true);
	// Pad entries for the mTypeRiemannSolver vector.
  PadEntriesVectorData(buffervec, "RIEMANN_SOLVER", mNZone, 1, 2);
	// Deduce the type from the buffer, based on the mapping.
	mTypeRiemannSolver = GenericVectorMap(MapTypeRiemannSolver, buffervec, "RIEMANN_SOLVER");

	// Read the mesh format specified.
	NInputUtility::AddScalarOption(paramfile, "MESH_FORMAT", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mMeshFormat = GenericScalarMap(MapMeshFormat, buffer, "MESH_FORMAT");

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


	// Read the writing frequency of the visualization files.
	NInputUtility::AddScalarOption(paramfile, "WRITE_VIS_FREQ", mWriteVisFreq, true);
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
	// Read the temporal step specified.
	NInputUtility::AddScalarOption(paramfile, "TIME_STEP", mTimeStep, true);
	// Read the maximum number of temporal iterations.
	NInputUtility::AddScalarOption(paramfile, "MAX_ITER_TIME", mMaxIterTime, true);
	// Read the CFL number for temporal stability.
	NInputUtility::AddScalarOption(paramfile, "CFL_NUMBER", mCFL, true);

	// Ensure the CFL makes sense.
	if( mCFL < C_ZERO ) ERROR("CFL numbner must be positive.");

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
	
	// Check for and extract interface/periodic boundary condition information, if present.
	ExtractInfoInterfaceBC(filename);

	// Check that all the boundary + interface markers are indeed unique. 
	if( mMarkerTag.size() )
	{
		// Ensure both interface and boundary markers are correct in size.
		if( mMarkerTag.size() != (mBoundaryParamMarker.size() + 2*mInterfaceParamMarker.size()) )
		{
			ERROR("Markers are not consistent in size.");
		}

		// Ensure the markers are indeed unique.
		for(size_t i=0; i<mMarkerTag.size(); i++)
		{
			for(size_t j=i+1; j<mMarkerTag.size(); j++)
			{
				if( mMarkerTag[i].first == mMarkerTag[j].first ) 
				{
					ERROR("Markers are not unique."); 
				}
			}
		}
	}
	else
	{
		// No boundary information is detected, this is nonsense.
		ERROR("Marker boundary/interface conditions are missing.");
	}

	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
}

//-----------------------------------------------------------------------------------

void CConfig::ExtractInfoInterfaceBC
(
 const char *filename
)
 /*
	* Function that reads the interface/periodic BC information, if present.
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
		NInputUtility::AddVectorOption(paramfile, "MARKER_BC_INTERFACE", buffer, defval, false);

		// Check if any new values are found.
		if( buffer != defval )
		{
			// Ensure the values specified are 4: iName, jName, rx, ry.
			if( buffer.size() != 4 )
			{
				ERROR("Interface BCs must specify two markers and a 2D translation vector.");
			}

			// Accumulate the information.
			mInterfaceParamMarker.push_back( std::make_unique<CInterfaceParamMarker>(buffer) );
		}
		else
		{
			// Prepare to exit from this loop.
			condition = false;
		}
	}

	// Accumulate the marker tag information, if need be.
	if( mInterfaceParamMarker.size() )
	{
		for( auto& interface: mInterfaceParamMarker )
		{
			// Accumulate both pair of boundaries.
			mMarkerTag.push_back( {interface->mName,         interface->GetTypeBC()} );
			mMarkerTag.push_back( {interface->mNameMatching, interface->GetTypeBC()} );
		}
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





