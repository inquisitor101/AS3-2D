#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometry;
class CMarker;


struct CInternalFaceGeometry
{
	size_t mIndexElementM;
	size_t mIndexElementP;
	
	unsigned short mIndexZone;
};

struct CBoundaryFaceGeometry
{
	size_t mIndexElement;
	size_t mIndexZone;

	EFaceLocation mFaceLocation;
};

struct CInterfaceFaceGeometry
{
	size_t mIndexElementM;
	size_t mIndexElementP;
	
	unsigned short mIndexZoneM; // TODO: remove these, since they are in CInterfaceGroup  
	unsigned short mIndexZoneP; 

	EFaceLocation mFaceLocationM; // TODO: remove these, since they are in CInterfaceGroup
	EFaceLocation mFaceLocationP;
};

struct CInterfaceGroup
{
  unsigned short mIndexZoneM;
  unsigned short mIndexZoneP;

  EFaceLocation mFaceLocationM;
  EFaceLocation mFaceLocationP;

  std::string mFaceNameM;
  std::string mFaceNameP;
  
  size_t mIndexBegin;
  size_t mIndexEnd;

  size_t GetnElem(void) const { return mIndexEnd - mIndexBegin; }
  size_t GetFaceIndex(size_t i) const { return mIndexBegin + i; }
};


class CMultizoneFaceGeometry
{
	public:
    CMultizoneFaceGeometry(ETypeDirection direction) : mDirection(direction) {}

		void InitializeInternalFaces(const CGeometry *geometry_container);

		void InitializeInterfaceFaces(const CConfig   *config_container,
				                          const CGeometry *geometry_structure);


		ETypeFaceGeometry GetFaceTypeFromIndex(size_t i) const
		{
			if( i < GetnInternalFaces() )                       return ETypeFaceGeometry::INTERNAL;
			if( i < GetnInternalFaces() + GetnBoundaryFaces() ) return ETypeFaceGeometry::BOUNDARY;
			if( i < GetnFacesTotal() )                          return ETypeFaceGeometry::INTERFACE;

			ERROR("Invalid face index specified.");
		}

		const as3vector1d<CInternalFaceGeometry>  &GetInternalFaces(void)  const { return mInternalFaces; }
		const as3vector1d<CBoundaryFaceGeometry>  &GetBoundaryFaces(void)  const { return mBoundaryFaces; }
		const as3vector1d<CInterfaceFaceGeometry> &GetInterfaceFaces(void) const { return mInterfaceFaces; }

		CInternalFaceGeometry& GetInternalFace(size_t i)
		{
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInternalFaces[i]; 
		}
		const CInternalFaceGeometry& GetInternalFace(size_t i) const 
		{
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInternalFaces[i]; 
		}

		CBoundaryFaceGeometry& GetBoundaryFace(size_t i)
		{ 
#if DEBUG
			if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mBoundaryFaces[i]; 
		}
		const CBoundaryFaceGeometry& GetBoundaryFace(size_t i) const 
		{ 
#if DEBUG
			if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mBoundaryFaces[i]; 
		}

		CInterfaceFaceGeometry& GetInterfaceFace(size_t i)
		{
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInterfaceFaces[i]; 
		}
		const CInterfaceFaceGeometry& GetInterfaceFace(size_t i) const 
		{
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
			return mInterfaceFaces[i]; 
		}

		size_t GetnInternalFaces(void)  const { return mInternalFaces.size();  }
		size_t GetnBoundaryFaces(void)  const { return mBoundaryFaces.size();  }
		size_t GetnInterfaceFaces(void) const { return mInterfaceFaces.size(); }

		size_t GetnFacesTotal(void) const
		{
			return GetnInternalFaces() + GetnBoundaryFaces() + GetnInterfaceFaces();
		}

    size_t GetnInterfaceGroups(void) const
    {
      return mInterfaceGroups.size();
    }

	private:
    ETypeDirection mDirection;

    as3vector1d<CInternalFaceGeometry>  mInternalFaces;
		as3vector1d<CBoundaryFaceGeometry>  mBoundaryFaces;
		as3vector1d<CInterfaceFaceGeometry> mInterfaceFaces;

    as3vector1d<CInterfaceGroup> mInterfaceGroups; // Family of interfaces, each containing the element faces on it.


    const CMarker* GetMatchingMarker(const CGeometry   *geometry_container,
                                     const std::string &marker_name);

    const EFaceLocation GetMarkerFaceLocation(const CMarker *marker_container);

		void CheckConformityMarkers(const CConfig         *config_container,
				                        const CGeometry       *geometry_container,
																const CMarker         *owner_marker,
																const CMarker         *match_marker,
																CInterfaceParamMarker *param_interface,
                                CInterfaceGroup       *interface_group);

    as3vector1d<size_t> GetIndexInterfacesAlongDirection(const CConfig   *config_container,
                                                         const CGeometry *geometry_container);

    ETypeDirection GetDirectionFromFaceLocation(EFaceLocation location) const
    {
      if( location == EFaceLocation::IMIN || location == EFaceLocation::IMAX ) return ETypeDirection::IDIR;
      if( location == EFaceLocation::JMIN || location == EFaceLocation::JMAX ) return ETypeDirection::JDIR;

      ERROR("Cannot deduce direction from face location.");
    }
};





