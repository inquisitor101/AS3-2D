#include "openmp_structure.hpp"


//-----------------------------------------------------------------------------------
// COpenMP member functions.
//-----------------------------------------------------------------------------------


COpenMP::COpenMP
(
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 as3vector1d<std::unique_ptr<ISolver>>    &solver_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container
)
 /*
	* Constructor for the OpenMP shared memory parallelization class.
	*/
{
	// Initialize the volume indices.
	InitializeVolumeIndices(config_container, geometry_container);

	// Initialize the internal face indices in the i-direction.
	InitializeSurfaceIFaces(config_container, geometry_container, interface_container);

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
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container
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



  // TESTING.
  size_t nIFaces = 0;
  for(auto& zone: geometry_container->GetZoneGeometry() )
  {
    size_t niElem = zone->GetnxElem();
    size_t njElem = zone->GetnyElem();
    nIFaces += njElem * (niElem+1);
  }

  mFacesIDir.resize(nIFaces);

  size_t IDX = 0;
  for(unsigned short iZone=0; iZone<nZone; iZone++)
  {
    const unsigned int nxElem = geometry_container->GetZoneGeometry(iZone)->GetnxElem();
    const unsigned int nyElem = geometry_container->GetZoneGeometry(iZone)->GetnyElem();

    for(size_t j=0; j<nyElem; j++)
    {
      for(size_t i=0; i<nxElem+1; i++)
      {
        const size_t iFace = j*(nxElem+1) + i;
        mFacesIDir[IDX++] = std::make_unique<CIndexElement>(iZone, iFace);
      }
    }
  }







  // CHANGED: below is the new version which should consider only unique faces.



  // Unique faces only.
  size_t nIFacesTotal = 0;
  for(auto& zone : geometry_container->GetZoneGeometry() )
  {
    size_t niElem = zone->GetnxElem();
    size_t njElem = zone->GetnyElem();
    nIFacesTotal += njElem * (niElem+1);
  }

  // Subtract all the interface faces, since they contain two faces and we are
  // only interested in the unique faces (i.e. any one of the two would suffice).
  size_t nInterfaceElements = 0;
  for( auto& interface : interface_container )
  {
    if( (interface->GetjFace() == EFaceElement::IMIN) or (interface->GetjFace() == EFaceElement::IMAX) )
    {
      nInterfaceElements += interface->GetnElem();
    }
  }

  // Deduce the number of unique faces.
  const size_t nUniqueIFaces = nIFacesTotal - nInterfaceElements;

  if( nUniqueIFaces < 0 ) 
  {
    ERROR("Detected a non-positive number of unique faces in the IDir: " + std::to_string(nUniqueIFaces));
  }


  // Allocate the correct amount of memory.
  as3vector1d<CFaceIndices> mUniqueIFaces;
  mUniqueIFaces.reserve(nUniqueIFaces);


  // Create a lambda function that detects whether a given face in the i-direction is unique or not.
  auto lGetisUniqueFaceIDir = [&](const size_t iZone,
                                  const size_t iFace,
                                  const size_t jFace,
                                  const EFaceElement face_pos) -> bool
  {
    const size_t iZone  = zone->GetZoneID();
    const size_t niElem = zone->GetnxElem();

    for (const auto& interface : interface_container)
    {
      // First, check if this face is part of an interface.
      // We only consider the jth face of this interface, since
      // the iFace is treated as the unique face.
      if (interface->GetjFace() != face_pos ||
          interface->GetjZone() != iZone)
      {
        continue;
      }
      
      for (const auto& [I, J] : interface->GetIndexElement() )
      {
        // Deduce the element contained on this interface's jth face.
        size_t iElemInterface = J % (niElem+1); 
        if(face_pos == EFaceElement::IMAX) iElemInterface--; // Correct for the IMAX face's element, since its on the left.

        // Check whether it corresponds to the same element.
        if (iElemInterface == J)
        {
          return false;
        }
      }
    }
    
    // If no matching interface/element was found, the face is unique.
    return true;
  };


  // Let's define the unique faces in the I-direction.
  size_t index = 0;
  for( const auto& zone : geometry_container->GetZoneGeometry() )
  {
    const size_t iZone  = zone->GetZoneID();
    const size_t niElem = zone->GetnxElem();
    const size_t njElem = zone->GetnyElem();

    for(auto jFace=0; jFace<njElem; jFace++)
    {
      // Internal faces are always included, since we treat them uniquely.
      for(auto iFace=1; iFace<niElem; iFace++)
      {
        const size_t iFaceFlattened = jFace * (niElem+1) + iFace;
        mUniqueIFaces.emplace_back( iZone, iFaceFlattened );

        if( lGetisUniqueFaceIDir(zone, iFace, jFace) )
        {
          const size_t iFaceFlattened = jFace * (niElem+1) + iFace;
          mUniqueIFaces.emplace_back( iZone, iFaceFlattened );
          index++;
        }
      }

      // Next are the IMIN faces, which might be a boundary or an interface.
      const size_t iFace = 0;
      {
        const size_t iFaceFlattened = jFace * (niElem+1) + iFace;
        if( lGetisUniqueFaceIDir(zone, iFace, jFace, EFaceElement::IMIN) )
        {
          const size_t iFaceFlattened = jFace * (niElem+1) + iFace;
          mUniqueIFaces.emplace_back( iZone, iFaceFlattened );
          index++;
        }
      }

      // TODO: IMAX face

    }
  }

  // Consistency check.
  if( mUniqueIFaces.size() != nUniqueIFaces )
  {
    std::cout << "mUniqueIFaces.size(): " << mUniqueIFaces.size() << ", nUniqueIFaces: " << nUniqueIFaces << std::endl;
    //ERROR("Encountered wrong number of unique i-faces.");
  }

  std::cout << "size of mUniqueIFaces: " << mUniqueIFaces.size() << ", index: " << index << ", nIFacesTotal: " << nIFacesTotal << ", nInterfaceElements: " << nInterfaceElements << std::endl;

  ERROR("");
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



