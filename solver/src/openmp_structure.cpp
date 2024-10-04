#include "openmp_structure.hpp"


//-----------------------------------------------------------------------------------
// COpenMP member functions.
//-----------------------------------------------------------------------------------


COpenMP::COpenMP
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
 /*
	* Constructor for the OpenMP shared memory parallelization class.
	*/
{
	// Initialize the volume indices.
	InitializeVolumeIndices(config_container, geometry_container);

	// Initialize the internal face indices in the i-direction.
	InitializeSurfaceIFaces(config_container, geometry_container);

	// Initialize the internal face indices in the j-direction.
	InitializeSurfaceJFaces(config_container, geometry_container);


	mResIMin.resize(mInternIFace.size());
	for(size_t i=0; i<mResIMin.size(); i++)
	{
		const unsigned short iZone  = mInternIFace[i]->mZone;
		const unsigned short nVar   = solver_container[iZone]->GetnVar();
		const unsigned short nSol2D = solver_container[iZone]->GetStandardElement()->GetnSol2D();
		mResIMin[i].resize(nVar, nSol2D);
	}


	mResJMin.resize(mInternJFace.size());
	for(size_t i=0; i<mResJMin.size(); i++)
	{
		const unsigned short iZone  = mInternJFace[i]->mZone;
		const unsigned short nVar   = solver_container[iZone]->GetnVar();
		const unsigned short nSol2D = solver_container[iZone]->GetStandardElement()->GetnSol2D();
		mResJMin[i].resize(nVar, nSol2D);
	}
}

//-----------------------------------------------------------------------------------

COpenMP::~COpenMP
(
 void
)
 /*
	* Destructor, which cleans up after the OpenMP class.
	*/
{

}

//-----------------------------------------------------------------------------------

void COpenMP::InitializeVolumeIndices
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the volume indices for each element. 
	*/
{
	// Get the number of zones.
	const unsigned short nZone = config_container->GetnZone();

	// Deduce the overall number of elements in all zones.
	size_t nElem = 0;
	for( auto& zone: geometry_container->GetZoneGeometry() )
	{
		nElem += zone->GetnElem();
	}

	// Allocate the memory for the volume indices.
	mIndexVolume.resize(nElem);

	// Map the indices, according to the loop convention used.
	size_t II = 0;
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		const unsigned int nElem = geometry_container->GetZoneGeometry(iZone)->GetnElem();
		for(unsigned int iElem=0; iElem<nElem; iElem++)
		{
			mIndexVolume[II++] = std::make_unique<CIndexElement>(iZone, iElem);
		}
	}
}

//-----------------------------------------------------------------------------------

void COpenMP::InitializeSurfaceIFaces
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the internal face indices in the i-direction. 
	*/
{
	// Get the number of zones.
	const unsigned short nZone = config_container->GetnZone();

	// Deduce the overall number of elements in all zones.
	size_t nElem = 0;
	for( auto& zone: geometry_container->GetZoneGeometry() )
	{
		nElem += zone->GetnElem() - zone->GetnyElem();
	}

	// Allocate the memory for the (right) i-surface indices.
	mInternIFace.resize(nElem);


	// Map the indices, according to the loop convention used.
	size_t II = 0;
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		const unsigned int nxElem = geometry_container->GetZoneGeometry(iZone)->GetnxElem();
		const unsigned int nyElem = geometry_container->GetZoneGeometry(iZone)->GetnyElem();
		
		for(size_t jElem=0; jElem<nyElem; jElem++)
		{
			for(size_t iElem=1; iElem<nxElem; iElem++)
			{
				// Right-element index (IMAX).
				const unsigned int iElemIMAX = jElem*nxElem + iElem;
				mInternIFace[II++] = std::make_unique<CIndexElement>(iZone, iElemIMAX);
			}
		}
	}
}

//-----------------------------------------------------------------------------------

void COpenMP::InitializeSurfaceJFaces
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the internal face indices in the j-direction. 
	*/
{
	// Get the number of zones.
	const unsigned short nZone = config_container->GetnZone();

	// Deduce the overall number of elements in all zones.
	size_t nElem = 0;
	for( auto& zone: geometry_container->GetZoneGeometry() )
	{
		nElem += zone->GetnElem() - zone->GetnyElem();
	}

	// Allocate the memory for the internal (top) j-surface indices.
	mInternJFace.resize(nElem);

	// Map the indices, according to the loop convention used.
	size_t II = 0;
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		const unsigned int nxElem = geometry_container->GetZoneGeometry(iZone)->GetnxElem();
		const unsigned int nyElem = geometry_container->GetZoneGeometry(iZone)->GetnyElem();
		
		for(size_t jElem=1; jElem<nyElem; jElem++)
		{
			for(size_t iElem=0; iElem<nxElem; iElem++)
			{
				// Top-element index (JMAX).
				const unsigned int iElemJMAX = jElem*nxElem + iElem;
				mInternJFace[II++] = std::make_unique<CIndexElement>(iZone, iElemJMAX);
			}
		}
	}
}



