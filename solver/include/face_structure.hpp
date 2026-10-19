#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "marker_structure.hpp"


// Forward declaration to avoid compiler problems.
class CMultizoneGeometry;
class CSinglezoneGeometry;




struct CElementFaceIndex
{
  std::size_t mIndexFamily;
  std::size_t mIndexFace;
};




//-----------------------------------------------------------------------------------


struct CInternalElementFaceGeometry
{
  std::size_t mIndexElementM;
};

//-----------------------------------------------------------------------------------

class CInternalFacesFamily
{
  using AInternalFaceVector = as3vector1d<CInternalElementFaceGeometry>;
  
  public:
    CInternalFacesFamily(ETypeFace                  face_type, 
                         const CSinglezoneGeometry *singlezone_geometry_container);

    const AInternalFaceVector& GetInternalFaces(void) const
    {
      return mInternalFaces;
    }
    const CInternalElementFaceGeometry& GetInternalFace(std::size_t i) const
    {
      return mInternalFaces[i];
    }
    std::size_t GetnFace(void) const {return mInternalFaces.size();}
    std::size_t GetiZone(void) const {return mIndexZone;}

  private:
    ETypeFace           mTypeFace;
    std::size_t         mIndexZone;
    std::size_t         mOffsetElementPlus;
    AInternalFaceVector mInternalFaces;

    void InitializeInternalFaces(const CSinglezoneGeometry *singlezone_geometry_container);
};

//-----------------------------------------------------------------------------------

struct CBoundaryElementFaceGeometry
{
  std::size_t mIndexElement;
};

//-----------------------------------------------------------------------------------

class CBoundaryFacesFamily
{
  using ABoundaryFaceVector = as3vector1d<CBoundaryElementFaceGeometry>;
  
  public:
    CBoundaryFacesFamily(const CSinglezoneGeometry                    *singlezone_geometry_container,
                         const CExternalFamilyMarker::CExternalMarker &marker,
                         const std::string                            &name);

    const ABoundaryFaceVector& GetBoundaryFaces(void) const
    {
      return mBoundaryFaces;
    }
    const CBoundaryElementFaceGeometry& GetBoundaryFace(std::size_t i) const
    {
      return mBoundaryFaces[i];
    }
    std::size_t GetnFace(void) const {return mBoundaryFaces.size();}
    std::size_t GetiZone(void) const {return mIndexZone;}
    std::string GetName(void) const {return mName;}
    EFaceLocation GetFaceLocation(void) const {return mFaceLocation;}
    bool GetisMaxFace(void) const {return mIsMaxFace;}

  private:
    ETypeFace           mTypeFace;
    EFaceLocation       mFaceLocation;
    std::size_t         mIndexZone;
    std::string         mName;
    bool                mIsMaxFace;
    ABoundaryFaceVector mBoundaryFaces;

    void InitializeBoundaryFaces(const CSinglezoneGeometry *singlezone_geometry_container);
};

//-----------------------------------------------------------------------------------

struct CInterfaceElementFaceGeometry
{
  std::size_t mIndexElementI;
  std::size_t mIndexElementJ;
};

//-----------------------------------------------------------------------------------

class CInterfaceFacesFamily
{
  using AInterfaceFaceVector = as3vector1d<CInterfaceElementFaceGeometry>;
  
  public:
    CInterfaceFacesFamily(const CMultizoneGeometry                     *multizone_geometry_container,
                          const CInternalFamilyMarker::CInternalMarker  &marker);

    CInterfaceFacesFamily(const CMultizoneGeometry                     *multizone_geometry_container,
                          const CPeriodicFamilyMarker::CPeriodicMarker  &marker);

    const AInterfaceFaceVector& GetInterfaceFaces(void) const
    {
      return mInterfaceFaces;
    }
    const CInterfaceElementFaceGeometry& GetInterfaceFace(std::size_t i) const
    {
      return mInterfaceFaces[i];
    }
    EFaceLocation GetFaceLocationFaceI(void) const {return mFaceLocationI;}
    std::size_t GetnFace(void) const {return mInterfaceFaces.size();}

    std::size_t GetiZone(void) const {return mIndexZoneI;}
    std::size_t GetjZone(void) const {return mIndexZoneJ;}

    EFaceLocation GetiFaceLocation(void) const {return mFaceLocationI;}
    EFaceLocation GetjFaceLocation(void) const {return mFaceLocationJ;}

		bool GetisMaxFaceI(void) const {return mIsMaxFaceI;}
		bool GetisMaxFaceJ(void) const {return mIsMaxFaceJ;}

    const std::string& GetiFaceName() const { return mFaceNameI; }
    const std::string& GetjFaceName() const { return mFaceNameJ; }

    ETypeFace GetiTypeFace(void) const { return mTypeFaceI; }
    ETypeFace GetjTypeFace(void) const { return mTypeFaceJ; }

  private:
    std::size_t mIndexZoneI;
    std::size_t mIndexZoneJ;

    ETypeFace mTypeFaceI;
    ETypeFace mTypeFaceJ;

    EFaceLocation mFaceLocationI;
    EFaceLocation mFaceLocationJ;

    std::string mFaceNameI;
    std::string mFaceNameJ;

		bool mIsMaxFaceI;
		bool mIsMaxFaceJ;

    bool mIsReversed; // should the j-edge be reversed?

    bool mIsPeriodic;

    std::array<as3double, 2> mTranslationVector;

    AInterfaceFaceVector mInterfaceFaces;

    void InitializeInterfaceFaces(const CMultizoneGeometry *multizone_geometry_container);
		void CheckConformityMarkers(const CMultizoneGeometry *multizone_geometry_container);
};

//-----------------------------------------------------------------------------------


// TODO: change this to CGroupFamilyFaces ?
template<typename TFamily>
class CGroupFaces
{
  private:
  
    struct CFaceRange
    {
      std::size_t GetnFace(void)   const { return mEnd - mBegin; }
      bool Contains(std::size_t i) const { return i >= mBegin && i < mEnd; }
      
      std::size_t mBegin;
      std::size_t mEnd;
    };
  
  public:
  
    CGroupFaces(void) = default;
  
    void Initialize(as3vector1d<TFamily>&& families)
    {
      mFamilies = std::move(families);

      mFamilyRanges.clear();
      mFamilyRanges.reserve( mFamilies.size() );

      std::size_t iGlobal = 0;    
      for(const auto& family : mFamilies)
      {
        const std::size_t nFace = family.GetnFace();
        mFamilyRanges.emplace_back( CFaceRange{iGlobal, iGlobal+nFace} );
        iGlobal += nFace;
      }
      
      mBegin = 0;
      mEnd   = iGlobal;
    }

    CElementFaceIndex FindElementFaceIndex(std::size_t iGlobal) const
    {
#if DEBUG
      if( iGlobal < mBegin || iGlobal >= mEnd ) ERROR("Invalid global face index.");
#endif
      const std::size_t iFamily    = FindIndexFamily(iGlobal);
      const CFaceRange& range = mFamilyRanges[iFamily];
  
      return CElementFaceIndex{ iFamily, iGlobal - range.mBegin };
    }
  
    std::size_t GetnFace(void) const
    {
#if DEBUG
      if( mBegin > mEnd ) ERROR("Invalid number of faces detected.");
#endif
      return mEnd - mBegin;
    }

    std::size_t GetnFamily(void) const
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

    const TFamily& GetFamily(std::size_t iFamily) const
    {
#if DEBUG
      if( iFamily >= mFamilies.size() ) ERROR("Invalid family index.");
#endif
      return mFamilies[iFamily];
    }
    TFamily& GetFamily(std::size_t iFamily)
    {
#if DEBUG
      if( iFamily >= mFamilies.size() ) ERROR("Invalid family index.");
#endif
      return mFamilies[iFamily];
    }
  
  private:
    std::size_t mBegin = 0;
    std::size_t mEnd   = 0;
 
    as3vector1d<TFamily>    mFamilies;
    as3vector1d<CFaceRange> mFamilyRanges;
  
    std::size_t FindIndexFamily(std::size_t iGlobal) const
    {
      for(std::size_t iFamily=0; iFamily<mFamilyRanges.size(); iFamily++)
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

		void InitializeFaces(const CConfig            *config_container,
				                 const CMultizoneGeometry *multizone_geometry_container);


		ETypeFaceGeometry DeduceFaceTypeFromIndex(std::size_t i) const
		{
			if( i < GetnInternalFaces() )                       return ETypeFaceGeometry::INTERNAL;
			if( i < GetnInternalFaces() + GetnBoundaryFaces() ) return ETypeFaceGeometry::BOUNDARY;
			if( i < GetnFacesTotal() )                          return ETypeFaceGeometry::INTERFACE;
			ERROR("Invalid face index specified.");
		}

    std::size_t GetIndexInternalFace(std::size_t iGlobal) const
    {
      return iGlobal;
    }
    std::size_t GetIndexBoundaryFace(std::size_t iGlobal) const
    {
      return iGlobal - GetnInternalFaces();
    }
    std::size_t GetIndexInterfaceFace(std::size_t iGlobal) const
    {
      return iGlobal - GetnInternalFaces() - GetnBoundaryFaces();
    }

    const auto& GetInternalFacesGroup(void)  const {return mInternalFacesGroup;}
    const auto& GetBoundaryFacesGroup(void)  const {return mBoundaryFacesGroup;}
    const auto& GetInterfaceFacesGroup(void) const {return mInterfaceFacesGroup;}

		std::size_t GetnInternalFaces(void)  const { return mInternalFacesGroup.GetnFace();  }
		std::size_t GetnBoundaryFaces(void)  const { return mBoundaryFacesGroup.GetnFace();  }
		std::size_t GetnInterfaceFaces(void) const { return mInterfaceFacesGroup.GetnFace(); }

		std::size_t GetnFacesTotal(void) const
		{
			return GetnInternalFaces() + GetnBoundaryFaces() + GetnInterfaceFaces();
		}

    ETypeFace GetTypeFace(void) const { return mTypeFace; }

    std::size_t GetnInterfaceFamilies(void) const { return mInterfaceFacesGroup.GetnFamily(); }

    // NOTE, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value! 
    CElementFaceIndex FindInternalElementFaceIndex(std::size_t i) const
    {
#if DEBUG
			if( i >= GetnInternalFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInternalFacesGroup.FindElementFaceIndex(i);
    } 

    const CInternalFacesFamily& GetInternalFacesFamily(std::size_t i) const
    {
#if DEBUG
			if( i >= mInternalFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInternalFacesGroup.GetFamily(i);
    }

    // NOTE, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value! 
    CElementFaceIndex FindInterfaceElementFaceIndex(std::size_t i) const
    {
#if DEBUG
			if( i >= GetnInterfaceFaces() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInterfaceFacesGroup.FindElementFaceIndex(i);
    } 
    
    const CInterfaceFacesFamily& GetInterfaceFacesFamily(std::size_t i) const
    {
#if DEBUG
			if( i >= mInterfaceFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size."); 
#endif
      return mInterfaceFacesGroup.GetFamily(i);
    }

    // Note, we return by value because its only 2 size_t variables, besides, the 
    // construction in CGroupFaces returns it by value!
    CElementFaceIndex FindBoundaryElementFaceIndex(std::size_t i) const
    {
#if DEBUG
      if( i >= GetnBoundaryFaces() ) ERROR("Index exceeds maximum data size.");
#endif
      return mBoundaryFacesGroup.FindElementFaceIndex(i);
    }

    const CBoundaryFacesFamily& GetBoundaryFacesFamily(std::size_t i) const
    {
#if DEBUG
      if( i >= mBoundaryFacesGroup.GetnFamily() ) ERROR("Index exceeds maximum data size.");
#endif
      return mBoundaryFacesGroup.GetFamily(i);
    }

	private:
    ETypeFace mTypeFace;

		// A Group is made up of: [iFamily][iFace].
    CGroupFaces<CInternalFacesFamily>  mInternalFacesGroup;
    CGroupFaces<CBoundaryFacesFamily>  mBoundaryFacesGroup;
    CGroupFaces<CInterfaceFacesFamily> mInterfaceFacesGroup;


		void InitializeInternalFaces(const CMultizoneGeometry *multizone_geometry_container);
		void InitializeBoundaryFaces(const CConfig            *config_container,
                                 const CMultizoneGeometry *multizone_geometry_container);
		void InitializeInterfaceFaces(const CConfig            *config_container,
				                          const CMultizoneGeometry *multizone_geometry_container);
};





