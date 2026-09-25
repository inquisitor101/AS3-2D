#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometry;
class CZoneGeometry;
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








struct CInternalElementFaceGeometry
{
  size_t mIndexElementM;
};

class CInternalFacesFamily
{
  using AInternalFaceVector = as3vector1d<CInternalElementFaceGeometry>;
  
  public:
    CInternalFacesFamily(ETypeFace            face_type, 
                         const CZoneGeometry *zone_geometry);

    const AInternalFaceVector& GetInternalFaces(void) const
    {
      return mInternalFaces;
    }
    const CInternalElementFaceGeometry& GetInternalFace(size_t i) const
    {
      return mInternalFaces[i];
    }
    size_t GetnFace(void) const {return mInternalFaces.size();}
    unsigned short GetiZone(void) const {return mIndexZone;}

  private:
    ETypeFace           mTypeFace;
    unsigned short      mIndexZone;
    size_t              mOffsetElementPlus;
    AInternalFaceVector mInternalFaces;

    void InitializeInternalFaces(const CZoneGeometry *zone_geometry);
};


struct CBoundaryElementFaceGeometry
{
  size_t mIndexElement;
};

class CBoundaryFacesFamily
{
  using ABoundaryFaceVector = as3vector1d<CBoundaryElementFaceGeometry>;
  
  public:
    CBoundaryFacesFamily(ETypeFace            face_type,
                         EFaceLocation        face_location,
                         const CZoneGeometry *zone_geometry);

    const ABoundaryFaceVector& GetBoundaryFaces(void) const
    {
      return mBoundaryFaces;
    }
    const CBoundaryElementFaceGeometry& GetBoundaryFace(size_t i) const
    {
      return mBoundaryFaces[i];
    }
    size_t GetnFace(void) const {return mBoundaryFaces.size();}
  
  private:
    ETypeFace           mTypeFace;
    EFaceLocation       mFaceLocation;
    unsigned short      mIndexZone;
    ABoundaryFaceVector mBoundaryFaces;
};


struct CInterfaceElementFaceGeometry
{
  size_t mIndexElementI;
  size_t mIndexElementJ;
};

class CInterfaceFacesFamily
{
  using AInterfaceFaceVector = as3vector1d<CInterfaceElementFaceGeometry>;
  
  public:
    CInterfaceFacesFamily(const CGeometry       *geometry_container,
                          const CMarker         *imarker_container,
                          const CMarker         *jmarker_container,
                          CInterfaceParamMarker *param_interface);

    void InitializeInterfaceFaces(const CGeometry       *geometry_container,
                                  const CMarker         *imarker_container,
                                  const CMarker         *jmarker_container,
                                  CInterfaceParamMarker *param_interface);

    const AInterfaceFaceVector& GetInterfaceFaces(void) const
    {
      return mInterfaceFaces;
    }
    const CInterfaceElementFaceGeometry& GetInterfaceFace(size_t i) const
    {
      return mInterfaceFaces[i];
    }
    EFaceLocation GetFaceLocationFaceI(void) const {return mFaceLocationI;}
    size_t GetnFace(void) const {return mInterfaceFaces.size();}

    unsigned short GetiZone(void) const {return mIndexZoneI;}
    unsigned short GetjZone(void) const {return mIndexZoneJ;}

    EFaceLocation GetiFaceLocation(void) const {return mFaceLocationI;}
    EFaceLocation GetjFaceLocation(void) const {return mFaceLocationJ;}

  private:
    ETypeFace   mTypeFaceI;
    ETypeFace   mTypeFaceJ;

    EFaceLocation mFaceLocationI;
    EFaceLocation mFaceLocationJ;

    std::string mFaceNameI;
    std::string mFaceNameJ;

    unsigned short mIndexZoneI;
    unsigned short mIndexZoneJ;

    AInterfaceFaceVector mInterfaceFaces;

		void CheckConformityMarkers(const CGeometry       *geometry_container,
																const CMarker         *imarker_container,
																const CMarker         *jmarker_container,
																CInterfaceParamMarker *param_interface);
};




struct CElementFaceIndex
{
  size_t mIndexFamily;
  size_t mIndexFace;
};

template<typename TFamily>
class CGroupFaces
{
  private:
  
    struct CFaceRange
    {
      size_t GetnFace(void)   const { return mEnd - mBegin; }
      bool Contains(size_t i) const { return i >= mBegin && i < mEnd; }
      
      size_t mBegin;
      size_t mEnd;
    };
  
  public:
  
    CGroupFaces(void) = default;
  
    void Initialize(as3vector1d<TFamily>&& families)
    {
      mFamilies = std::move(families);

      mFamilyRanges.clear();
      mFamilyRanges.reserve( mFamilies.size() );

      size_t iGlobal = 0;    
      for(const auto& family : mFamilies)
      {
        const size_t nFace = family.GetnFace();
        mFamilyRanges.emplace_back( CFaceRange{iGlobal, iGlobal+nFace} );
        iGlobal += nFace;
      }
      
      mBegin = 0;
      mEnd   = iGlobal;
    }

    CElementFaceIndex GetElementFaceIndex(size_t iGlobal) const
    {
#if DEBUG
      if( iGlobal < mBegin || iGlobal >= mEnd ) ERROR("Invalid global face index.");
#endif
      const size_t iFamily    = FindIndexFamily(iGlobal);
      const CFaceRange& range = mFamilyRanges[iFamily];
  
      return CElementFaceIndex{ iFamily, iGlobal - range.mBegin };
    }
  
    size_t GetnFace(void) const
    {
#if DEBUG
      if( mBegin > mEnd ) ERROR("Invalid number of faces detected.");
#endif
      return mEnd - mBegin;
    }

    size_t GetnFamily(void) const
    {
      return mFamilies.size();
    }

    const as3vector1d<TFamily>& GetFamilies(void) const
    {
      return mFamilies;
    }
    as3vector1d<TFamily>& GetFamilies(void)
    {
      return mFamilies;
    }

    const TFamily& GetFamily(size_t iFamily) const
    {
#if DEBUG
      if( iFamily >= mFamilies.size() ) ERROR("Invalid family index.");
#endif
      return mFamilies[iFamily];
    }
    TFamily& GetFamily(size_t iFamily)
    {
#if DEBUG
      if( iFamily >= mFamilies.size() ) ERROR("Invalid family index.");
#endif
      return mFamilies[iFamily];
    }
  
  private:
    size_t mBegin = 0;
    size_t mEnd   = 0;
 
    as3vector1d<TFamily>    mFamilies;
    as3vector1d<CFaceRange> mFamilyRanges;
  
    size_t FindIndexFamily(size_t iGlobal) const
    {
      for(size_t iFamily=0; iFamily<mFamilyRanges.size(); iFamily++)
      {
        if( mFamilyRanges[iFamily].Contains(iGlobal)) return iFamily;
      }
      ERROR("Invalid global face index.");
    }
};




class CMultizoneFaceGeometry
{
  // NOTE, this class always stores the global faces according to:
  //       [internal][boundary][interface]. 
	public:
    
    CMultizoneFaceGeometry(ETypeFace face_type) : mTypeFace(face_type) {}

		void InitializeFaces(const CConfig   *config_container,
				                 const CGeometry *geometry_structure);


		ETypeFaceGeometry GetFaceTypeFromIndex(size_t i) const
		{
			if( i < GetnInternalFaces() )                       return ETypeFaceGeometry::INTERNAL;
			if( i < GetnInternalFaces() + GetnBoundaryFaces() ) return ETypeFaceGeometry::BOUNDARY;
			if( i < GetnFacesTotal() )                          return ETypeFaceGeometry::INTERFACE;
			ERROR("Invalid face index specified.");
		}

    size_t GetIndexInternalFace(size_t iGlobal) const
    {
      return iGlobal;
    }
    size_t GetIndexBoundaryFace(size_t iGlobal) const
    {
      return iGlobal - GetnInternalFaces();
    }
    size_t GetIndexInterfaceFace(size_t iGlobal) const
    {
      return iGlobal - GetnInternalFaces() - GetnBoundaryFaces();
    }

    const auto& GetInternalFacesGroup(void)  const {return mInternalFacesGroup;}
    const auto& GetBoundaryFacesGroup(void)  const {return mBoundaryFacesGroup;}
    const auto& GetInterfaceFacesGroup(void) const {return mInterfaceFacesGroup;}


    // NOTE, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value! 
    CElementFaceIndex GetInternalElementFaceIndex(size_t i) const
    {
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInternalFacesGroup.GetElementFaceIndex(i);
    } 

    const CInternalFacesFamily& GetInternalFacesFamily(size_t i) const
    {
#if DEBUG
			if( i >= mInternalFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInternalFacesGroup.GetFamily(i);
    }




//		CBoundaryFaceGeometry& GetBoundaryFace(size_t i)
//		{ 
//#if DEBUG
//			if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size."); 
//#endif
//			return mBoundaryFaces[i]; 
//		}
//		const CBoundaryFaceGeometry& GetBoundaryFace(size_t i) const 
//		{ 
//#if DEBUG
//			if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size."); 
//#endif
//			return mBoundaryFaces[i]; 
//		}


    // NOTE, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value! 
    CElementFaceIndex GetInterfaceElementFaceIndex(size_t i) const
    {
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInterfaceFacesGroup.GetElementFaceIndex(i);
    } 
    
    const CInterfaceFacesFamily& GetInterfaceFacesFamily(size_t i) const
    {
#if DEBUG
			if( i >= mInterfaceFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInterfaceFacesGroup.GetFamily(i);
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

		size_t GetnInternalFaces(void)  const { return mInternalFacesGroup.GetnFace();  }
		size_t GetnBoundaryFaces(void)  const { return mBoundaryFacesGroup.GetnFace();  }
		size_t GetnInterfaceFaces(void) const { return mInterfaceFacesGroup.GetnFace(); }


		size_t GetnFacesTotal(void) const
		{
			return GetnInternalFaces() + GetnBoundaryFaces() + GetnInterfaceFaces();
		}

    size_t GetnInterfaceGroups(void) const
    {
      return mInterfaceGroups.size();
    }

    ETypeFace GetTypeFace(void) const { return mTypeFace; }

	private:
    
    ETypeFace mTypeFace;

    as3vector1d<CInternalFaceGeometry>  mInternalFaces;
		as3vector1d<CBoundaryFaceGeometry>  mBoundaryFaces;
		as3vector1d<CInterfaceFaceGeometry> mInterfaceFaces;

    as3vector1d<CInterfaceGroup> mInterfaceGroups; // Family of interfaces, each containing the element faces on it.


    // TODO: remove the previous versions and their respective functions and use the below.
    CGroupFaces<CInternalFacesFamily>  mInternalFacesGroup;
    CGroupFaces<CBoundaryFacesFamily>  mBoundaryFacesGroup;
    CGroupFaces<CInterfaceFacesFamily> mInterfaceFacesGroup;



    // TODO: CHANGE this or put it somewhere elegantly? or even use it!
    size_t GetFaceStartIndexOffset(ETypeFaceGeometry type) const
    {
      switch(type)
      {
        case(ETypeFaceGeometry::INTERNAL ): return 0;
        case(ETypeFaceGeometry::BOUNDARY ): return GetnInternalFaces();
        case(ETypeFaceGeometry::INTERFACE): return GetnInternalFaces() + GetnBoundaryFaces();
      }
      ERROR("Invalid face type.");
    }


		void InitializeInternalFaces(const CGeometry *geometry_container);

		void InitializeInterfaceFaces(const CConfig   *config_container,
				                          const CGeometry *geometry_structure);


    // TODO: remove this.
    const CMarker* GetMatchingMarker(const CGeometry   *geometry_container,
                                     const std::string &marker_name);

    // TODO: remove this.
    const EFaceLocation GetMarkerFaceLocation(const CMarker *marker_container);

		void CheckConformityMarkers(const CConfig         *config_container,
				                        const CGeometry       *geometry_container,
																const CMarker         *owner_marker,
																const CMarker         *match_marker,
																CInterfaceParamMarker *param_interface,
                                CInterfaceGroup       *interface_group);

    // TODO: remove this.
    as3vector1d<size_t> GetIndexInterfacesAlongDirection(const CConfig   *config_container,
                                                         const CGeometry *geometry_container);

    // TODO: remove this.
    ETypeFace GetFaceTypeFromFaceLocation(EFaceLocation location) const
    {
      if( location == EFaceLocation::IMIN || location == EFaceLocation::IMAX ) return ETypeFace::IFACE;
      if( location == EFaceLocation::JMIN || location == EFaceLocation::JMAX ) return ETypeFace::JFACE;

      ERROR("Cannot deduce face type from face location.");
    }
};





