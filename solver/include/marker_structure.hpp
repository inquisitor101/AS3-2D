#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "input_structure.hpp"

// Forward declaration to avoid compiler problems.
class CMultizoneGeometry;
class CExternalFamilyMarker;



namespace NMarkerUtility
{
  ETypeFace GetTypeFace(EFaceLocation face_location);
}



// NOTE, markers are divided into three:
//   *) external, which are boundary markers.
//   *) internal, which are interface markers.
//   *) periodic, which.. are periodic markers.
// ... all are defined on each zone's entire edge.


// Marker for a family of periodic zone-boundaries.
class CPeriodicFamilyMarker
{
  public:
    
    // Marker for a periodic zone-boundary.
    struct CPeriodicMarker
    {
      size_t mIndexZoneI;
      size_t mIndexZoneJ;

      EFaceLocation mFaceLocationI;
      EFaceLocation mFaceLocationJ;

      std::string mNameI;
      std::string mNameJ;
      
      bool mIsReversed;

      std::array<as3double, 2> mTranslationVector; // From I to J.
    
      // Return the face direction on each side of this marker.
      ETypeFace GetiTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationI);
      }

      ETypeFace GetjTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationJ);
      }
    };

    CPeriodicFamilyMarker(const CMultizoneGeometry       *multizone_geometry_container,
                          const CExternalFamilyMarker    &iexternal_marker,
                          const CExternalFamilyMarker    &jexternal_marker,
                          const std::array<as3double, 2> &translation);

    std::size_t GetnMarkers(void) const {return mMarkers.size();}
    const auto& GetMarkers(void) const {return mMarkers;}
    auto& GetMarkers(void) {return mMarkers;}
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}
  
  private:
    
    as3vector1d<CPeriodicMarker> mMarkers;
};



// Marker for a family of internal zone-boundaries (i.e. interfaces).
class CInternalFamilyMarker
{
  public:
    
    // Marker for an internal zone-boundary (i.e. interface).
    struct CInternalMarker
    {
      size_t mIndexZoneI;
      size_t mIndexZoneJ;

      EFaceLocation mFaceLocationI;
      EFaceLocation mFaceLocationJ;

      bool mIsReversed;

      // Return the face direction on each side of this marker.
      ETypeFace GetiTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationI);
      }

      ETypeFace GetjTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationJ);
      }
    };

    CInternalFamilyMarker(std::size_t nMarker)
    {
      mMarkers.reserve( nMarker );
    }

    std::size_t GetnMarkers(void) const {return mMarkers.size();}
    const auto& GetMarkers(void) const {return mMarkers;}
    auto& GetMarkers(void) {return mMarkers;}
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}
  
  private:
    
    as3vector1d<CInternalMarker> mMarkers;
};


// Marker for a family of external zone-boundaries.
class CExternalFamilyMarker
{
  public:

    // Marker for an external zone-boundary.
    struct CExternalMarker
    {
      std::size_t   mIndexZone;
      EFaceLocation mFaceLocation;

      // Return the face direction on this marker.
      ETypeFace GetTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocation);
      }
    };

    CExternalFamilyMarker(std::string name,
                          std::size_t nMarker)
      : mName(std::move(name))
    {
      mMarkers.reserve( nMarker );
    }

    std::size_t GetnMarkers(void) const {return mMarkers.size();}
    const std::string& GetName(void) const {return mName;}
    const auto& GetMarkers(void) const {return mMarkers;}
    auto& GetMarkers(void) {return mMarkers;}
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}

  private:
    
    std::string mName;
    as3vector1d<CExternalMarker> mMarkers;
};




