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

	// Read the zone connectivity filename.
	NInputUtility::AddScalarOption(paramfile, "ZONE_CONNECTIVITY_FILENAME", mZoneConnFilename, true);

	// Read the format of the input grid files.
	NInputUtility::AddScalarOption(paramfile, "INPUT_GRID_FORMAT", buffer, true);
	// Deduce the type from the buffer, based on the mapping.
	mInputGridFormat = GenericScalarMap(MapFormatFile, buffer, "INPUT_GRID_FORMAT");

	// Check if the extension is provided to the connectivity file, otherwise add it.
	if( ( NInputUtility::GetFileExtension(mZoneConnFilename) ).empty() )
	{
		mZoneConnFilename += ".cfg";
	}

	// Close file.
  paramfile.close();

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

	// Check what else to extract, based on the type of IC.
	switch(mTypeIC)
	{
		case(ETypeIC::GAUSSIAN_PRESSURE): { IC_GaussianPressure(filename); break; }
		case(ETypeIC::ISENTROPIC_VORTEX): { IC_IsentropicVortex(filename); break; }
		default: ERROR("Cannot extract additional information from the specified IC.");
	}


	// Close file.
  paramfile.close();

	// Report output.
	std::cout << " Done." << std::endl;

  // Return happily.
	return true;
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





