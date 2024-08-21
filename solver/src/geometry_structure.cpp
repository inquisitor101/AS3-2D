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
	// Check and report output for the existance of the specified grid files.
	CheckExistanceGridFiles(config_container);

	// Reserve memory for the grid zones.
	mZoneGeometry.resize(mNZone);

	// Indentify the information in each grid zone, separately.
	for(unsigned short iZone=0; iZone<mNZone; iZone++)
	{
		mZoneGeometry[iZone] = std::make_unique<CZoneGeometry>(config_container,
																						               config_container->GetZoneGridFilename(iZone),
																						               iZone);	
	}
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

void CGeometry::CheckExistanceGridFiles
(
 CConfig *config_container
)
 /*
	* Function that checks the existance of the grid files and reports the output.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n"
						<< "Generating geometry: " << std::endl;

	// Report output.
	std::cout << "  expecting " << mNZone << " grid files. " << std::endl;
	std::cout << "   ... detected: " << std::endl;

	// Loop over all the expected zones and report their grid filenames.
	for(unsigned short i=0; i<mNZone; i++)
	{
		std::string filename = config_container->GetZoneGridFilename(i);
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

	// Report output.
	std::cout << "Done." << std::endl;
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
 CConfig                   *config_container,
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

	// Allocate the correct number of markers in this zone.
	mMarker.resize( nMark );

	// Loop over the markers and initialize them.
	for(size_t i=0; i<nMark; i++)
	{
		// Determine the boundary condition associated with this marker. 
		ETypeBC bc = config_container->DetermineMarkerBC(name[i]); 

		// Initialize the marker.
		mMarker[i] = std::make_unique<CMarker>(mZoneID, name[i], bc, mark[i], face[i]); 
	}
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
