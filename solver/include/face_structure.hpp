#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometry;
class CZoneGeometry;
class CMarker;




struct CElementFaceIndex
{
  size_t mIndexFamily;
  size_t mIndexFace;
};




//-----------------------------------------------------------------------------------


struct CInternalElementFaceGeometry
{
  size_t mIndexElementM;
};

//-----------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------

struct CBoundaryElementFaceGeometry
{
  size_t mIndexElement;
};

//-----------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------

struct CInterfaceElementFaceGeometry
{
  size_t mIndexElementI;
  size_t mIndexElementJ;
};

//-----------------------------------------------------------------------------------

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

		bool GetisMaxFaceI(void) const {return mIsMaxFaceI;}
		bool GetisMaxFaceJ(void) const {return mIsMaxFaceJ;}

    const std::string& GetiFaceName() const { return mFaceNameI; }
    const std::string& GetjFaceName() const { return mFaceNameJ; }

    ETypeFace GetiTypeFace(void) const { return mTypeFaceI; }
    ETypeFace GetjTypeFace(void) const { return mTypeFaceJ; }

  private:
    ETypeFace   mTypeFaceI;
    ETypeFace   mTypeFaceJ;

    EFaceLocation mFaceLocationI;
    EFaceLocation mFaceLocationJ;

    std::string mFaceNameI;
    std::string mFaceNameJ;

    unsigned short mIndexZoneI;
    unsigned short mIndexZoneJ;

		bool mIsMaxFaceI;
		bool mIsMaxFaceJ;

    AInterfaceFaceVector mInterfaceFaces;

		void CheckConformityMarkers(const CGeometry       *geometry_container,
																const CMarker         *imarker_container,
																const CMarker         *jmarker_container,
																CInterfaceParamMarker *param_interface);
};

//-----------------------------------------------------------------------------------


// TODO: change this to CGroupFamilyFaces ?
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

    CElementFaceIndex FindElementFaceIndex(size_t iGlobal) const
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

//-----------------------------------------------------------------------------------


class CMultizoneFaceGeometry
{
  // NOTE, this class always stores the global faces according to:
  //       [internal][boundary][interface]. 
	public:
    
    CMultizoneFaceGeometry(ETypeFace face_type) : mTypeFace(face_type) {}

		void InitializeFaces(const CConfig   *config_container,
				                 const CGeometry *geometry_structure);


		ETypeFaceGeometry DeduceFaceTypeFromIndex(size_t i) const
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

		size_t GetnInternalFaces(void)  const { return mInternalFacesGroup.GetnFace();  }
		size_t GetnBoundaryFaces(void)  const { return mBoundaryFacesGroup.GetnFace();  }
		size_t GetnInterfaceFaces(void) const { return mInterfaceFacesGroup.GetnFace(); }

		size_t GetnFacesTotal(void) const
		{
			return GetnInternalFaces() + GetnBoundaryFaces() + GetnInterfaceFaces();
		}

    ETypeFace GetTypeFace(void) const { return mTypeFace; }

    size_t GetnInterfaceFamilies(void) const { return mInterfaceFacesGroup.GetnFamily(); }

    // NOTE, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value! 
    CElementFaceIndex FindInternalElementFaceIndex(size_t i) const
    {
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInternalFacesGroup.FindElementFaceIndex(i);
    } 

    const CInternalFacesFamily& GetInternalFacesFamily(size_t i) const
    {
#if DEBUG
			if( i >= mInternalFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInternalFacesGroup.GetFamily(i);
    }

    // NOTE, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value! 
    CElementFaceIndex FindInterfaceElementFaceIndex(size_t i) const
    {
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInterfaceFacesGroup.FindElementFaceIndex(i);
    } 
    
    const CInterfaceFacesFamily& GetInterfaceFacesFamily(size_t i) const
    {
#if DEBUG
			if( i >= mInterfaceFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInterfaceFacesGroup.GetFamily(i);
    }

	private:
    ETypeFace mTypeFace;

		// A Group is made up of: [iFamily][iFace].
    CGroupFaces<CInternalFacesFamily>  mInternalFacesGroup;
    CGroupFaces<CBoundaryFacesFamily>  mBoundaryFacesGroup;
    CGroupFaces<CInterfaceFacesFamily> mInterfaceFacesGroup;


		void InitializeInternalFaces(const CGeometry *geometry_container);

		void InitializeInterfaceFaces(const CConfig   *config_container,
				                          const CGeometry *geometry_structure);
};





