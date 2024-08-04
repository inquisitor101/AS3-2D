#include "geometry_structure.hpp"



//-----------------------------------------------------------------------------------
// CGeometry member functions.
//-----------------------------------------------------------------------------------


CGeometry::CGeometry
(
 CConfig *config_container
)
	:
		mNZone( config_container->GetnZone() )
 /*
	* Constructor for the geometry, which contains all the grid geometry.
	*/
{
	// Extract the information contained in the zone connectivity file.
	ExtractZoneConnectivity(config_container);

	// Reserve memory for the grid zones.
	mZoneGeometry.resize(mNZone);

	// Indentify the information in each grid zone, separately.
	for(unsigned short iZone=0; iZone<mNZone; iZone++)
	{
		mZoneGeometry[iZone] = std::make_unique<CZoneGeometry>(config_container,
																						               mZoneGridFilename[iZone],
																						               iZone);	
	}

	// Match each interface marker to its corresponding zone and tag.
	MatchInterfaceMarkers(config_container);
}

//-----------------------------------------------------------------------------------

CGeometry::~CGeometry
(
 void
)
 /*
	* Destructor, which cleans up after the driver class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CGeometry::ExtractZoneConnectivity
(
 CConfig *config_container
)
 /*
	* Function that extracts the grid zone connectivity information.
	*/
{
	// First, check if the connectivity has the right extension, otherwise add it.
	std::string filename = config_container->GetZoneConnFilename();

 	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Reading zone connectivitiy file: " << filename << std::endl;

	// Message stream.
	std::ostringstream message;

  // Check if file exists.
  std::ifstream connfile(filename);
  if( !connfile.good() ) ERROR("Connectivity file could not be opened!");

	// Report output.
	std::cout << "  expecting " << mNZone << " grid files. " << std::endl;
	std::cout << "  ... detected: " << std::endl;

	// Reserve memory for the grid filenames.
	mZoneGridFilename.resize(mNZone);

	// Loop over all the expected zones and read their grid filenames.
	for(unsigned short i=0; i<mNZone; i++)
	{
		// Current zone filename.
		std::string ifile = "GRID_FILENAME_ZONE_" + std::to_string(i);
		NInputUtility::AddScalarOption(connfile, ifile.c_str(), mZoneGridFilename[i], true);

		// Ensure that the grid extension is correct.
		switch( config_container->GetMeshFormat() )
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

		// Report output.
		std::cout << "      " << i << "): " << mZoneGridFilename[i] << "\n";
	}

	// Close file.
  connfile.close();

	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void CGeometry::MatchInterfaceMarkers
(
 CConfig *config_container
)
 /*
	* Function that maps each interface to its assigned neighbor.
	*/
{
	// Max length of marker names, used for output format.
	size_t maxlen = 0;

	// Loop over all the zones and match their data.
	for( auto& zone : mZoneGeometry )
	{
		// Extract the ineterface-type marker objects.
		auto& imark = zone->GetMarkerInterface();

		// Loop over each marker and assign its match.
		for( auto& im : imark )
		{
			// Get current zone index.
			auto  iZone = im->GetZoneID();

			// Get matching zone information.
			auto  jZone = im->GetMatchingZoneID();
			auto& jName = im->GetMatchingNameMarkerTag();

			// Ensure matching zone is different than current zone.
			if( iZone == jZone ) ERROR("Interface zones must be different.");

			// Update the max length of the marker name tags.
			maxlen = std::max( maxlen, jName.size() );

			// Extract the markers of the matching zone.
			auto& jmark = mZoneGeometry[jZone]->GetMarkerInterface();

			// Loop over the matching zone markers and search for this unique marker.
			bool Found = false;
			for( auto& jm : jmark )
			{
				// Get the current marker information.
				auto  zone = jm->GetZoneID();
				auto& name = jm->GetNameMarkerTag();
			
				// See if this marker is indeed the one we are interested in.
				if( (zone == jZone) && (name == jName)  )
				{
					// Match the current markers together.
					im->SetMatchingMarker(std::static_pointer_cast<CMarker>(jm));

					// Search is over, skip the loop.
					Found = true;
					break;
				}
			}

			// Check if the markers are found.
			if( !Found ) ERROR("Could not find the matching markers.");
		}
	}


 	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "The zone marker interfaces detected and paired are: " << std::endl;


	// Since all the markers so far are interface markers, they all should be shared once.
	for( auto& z : mZoneGeometry )
	{
		// Report the output.
		std::cout << "\n";

		for( auto& m : z->GetMarkerInterface() )
		{
			// Lock the weak pointer and check if its valid.
			auto ptr = m->GetMatchingMarker().lock();

			if( ptr )
			{
				// Copy the name tags to a fixed-length string.
				std::string tmp = m->GetNameMarkerTag();
				std::string tag(maxlen, ' ');
				for(size_t ii=0; ii<tmp.size(); ii++) tag[ii] = tmp[ii];
				std::cout << "  iZone: " << m->GetZoneID()
					        << ", iName: " << tag
									<< " ==>"
									<< "  jZone: " << ptr->GetZoneID()
									<< ", jName: " << ptr->GetNameMarkerTag()
									<< std::endl;
			}
			else
			{
				ERROR("Some markers could not be paired.");
			}
		}
	}

	// Report output.
	std::cout << "\nDone." << std::endl;
}


//-----------------------------------------------------------------------------------
// CZoneGeometry member functions.
//-----------------------------------------------------------------------------------


CZoneGeometry::CZoneGeometry
(
 CConfig        *config_container,
 std::string     gridfile,
 unsigned short  iZone
)
	: 
		mZoneID(iZone), 
	  mGridFile(gridfile),
		mNPolyGrid(config_container->GetnPolyGrid(iZone))
 /*
	* Constructor for the zone geometry, which contains a single grid geometry.
	*/
{
	// Read the marker interface information for this zone only.
	ReadMarkerInterfaceZone(config_container);
}

//-----------------------------------------------------------------------------------

CZoneGeometry::~CZoneGeometry
(
 void
)
 /*
	* Destructor, which cleans up after the driver class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CZoneGeometry::InitializeElements
(
 as3vector2d<double> &x,
 as3vector2d<double> &y,
 unsigned int         nxElem,
 unsigned int         nyElem
)
 /*
	* Function that initializes and defines all the element geometry in this zone..
	*/
{
	// Deduce the number of elements and ensure consistency.
	if( x.size() != y.size() ) ERROR("Number of coordinates in x and y is not identical.");

	// Ensure the total number of elements is correct.
	if( static_cast<size_t>(nxElem*nyElem) != x.size() ) ERROR("Inconsistency in number of elements.");

	// Set the number of elements in each dimension.
	mNxElem = nxElem;
	mNyElem = nyElem;

	// Allocate memory for the total elements.
	mElementGeometry.resize( x.size() );

	// Total number of grid nodes in 2D, as user-specified.
	const size_t nNode2D = (mNPolyGrid+1)*(mNPolyGrid+1);

	// Loop over each element and instantiate its coordinates.
	for(size_t i=0; i<x.size(); i++)
	{
		// Ensure the polynomial order is correct.
		if( (x[i].size() != nNode2D) || (y[i].size() != nNode2D) ) ERROR("Elements do not match polynomial order.");

		// Instantiate the current element.
		mElementGeometry[i] = std::make_unique<CElementGeometry>( x[i], y[i] );
	}
}

//-----------------------------------------------------------------------------------

void CZoneGeometry::InitializeMarkers
(
 as3vector2d<unsigned int> &mark,
 as3vector2d<EFaceElement> &face,
 as3vector1d<std::string>  &name
)
 /*
	* Function that defines all the interface markers in this zone..
	*/
{
	// For convenience, extract the number of total markers (interface+physical).
	const size_t nMark = mark.size();

	// Ensure the number of interface markers makes sense.
	if( ( mMarkerInterface.size() > nMark ) || ( name.size() != nMark ) ) ERROR("Inconsistent number of markers detected.");

	for( auto& im: mMarkerInterface )
	{
		// Flag, whether the marker is found or not.
		bool found = false;
	
		// Search for this marker and assign its properties.
		for(size_t j=0; j<nMark; j++)
		{
			auto& jn = name[j];
			if( im->GetNameMarkerTag() == jn )
			{
				// Ensure the marker is specified as an interface.
				if( im->GetTypeMarker() != ETypeBCs::INTERFACE ) ERROR("Incorrect marker type.");

				// Set its corresponding information.
				im->SetElementIndices(mark[j]);
				im->SetElementFaces(face[j]);

				found = true;
				break;
			}
		}
		if( !found ) ERROR("Could not locate the marker.");
	}
}

//-----------------------------------------------------------------------------------

void CZoneGeometry::ReadMarkerInterfaceZone
(
 CConfig *config_container
)
 /*
	* Function that reads the marker interface for this zone.
	*/
{
  // Open connectivity file.
  std::ifstream connfile( config_container->GetZoneConnFilename() );

  // Buffers for storing temporary strings.
	as3vector1d<std::string> buffer;
	unsigned int             nmark;

	// Temporary storage for the keyword.
	std::string keyword;

	
	// Specify number of interface markers in this zone.
	keyword = "NUMBER_INTERFACE_ZONE_" + std::to_string(mZoneID);
	// Read the number of markers.
	NInputUtility::AddScalarOption(connfile, keyword.c_str(), nmark, true);

	// Extract the interface marker tags in this zone.
	for(size_t iMarker=0; iMarker<nmark; iMarker++)
	{
		// Assign the current keyword to search for.
		keyword = "MARKER_INTERFACE_ZONE_" + std::to_string(mZoneID);

		// Read each marker specification.
		NInputUtility::AddVectorOption(connfile, keyword.c_str(), buffer, false);

		// Consistency check.
		if( buffer.size() != 3 ) ERROR("Interface marker information must contain 3 values.");

		// Individually extract the information, f.
		std::string nameI = buffer[0];
		std::string nameJ = buffer[2];

		// By convention, the matching zone index must be the last digit in the name (e.g. ZONE_#).
		unsigned short zoneJ = NInputUtility::ExtractLastDigitFromString<unsigned short>(buffer[1]);

		// Instantiate the associated marker object.
		mMarkerInterface.push_back( 
				std::make_shared<CInterfaceMarker>(mZoneID, nameI, zoneJ, nameJ, ETypeBCs::INTERFACE) );
	}


	// Close file.
  connfile.close();
}


//-----------------------------------------------------------------------------------
// CElementGeometry member functions.
//-----------------------------------------------------------------------------------


CElementGeometry::CElementGeometry
(
 as3vector1d<double> &x,
 as3vector1d<double> &y
)
 /*
	* Constructor for the element geometry, which contains a single element geometry.
	*/
{
	// Ensure consistency.
	if( x.size() != y.size() ) ERROR("Inconsistent size in coordinates.");

	// Allocate memory for the coordinates.
	mCoordSolDOFs.resize( 2, x.size() );

	// Copy the coordinates.
	for(size_t i=0; i<x.size(); i++)
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
	* Destructor, which cleans up after the element geometry class.
	*/
{

}
