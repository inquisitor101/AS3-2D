#pragma once

/*!
 * @file
 * @brief Zone-edge markers for physical, internal and periodic boundaries.
 */


#include <array>
#include <cstddef>
#include <string>
#include <utility>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "input_structure.hpp"

// Forward declaration to avoid compiler problems.
class CMultizoneGeometry;
class CExternalFamilyMarker;



/*!
 * @brief Conversion between zone-edge locations and logical face directions.
 */
namespace NMarkerUtility
{

  /*!
   * @brief Getter function which returns the logical face direction associated with the local edge.
   *
   * @param[in] face_location local edge: IMIN, IMAX, JMIN or JMAX.
   *
   * @return logical face direction associated with the local edge.
   */
  ETypeFace GetTypeFace(EFaceLocation face_location);
}



// NOTE, markers are divided into three:
//   *) external, which are boundary markers.
//   *) internal, which are interface markers.
//   *) periodic, which.. are periodic markers.
// ... all are defined on each zone's entire edge.


/*!
 * @brief Pairs of zone-edge markers connected by a periodic translation.
 */
class CPeriodicFamilyMarker
{
  public:
    

    /*!
     * @brief Two zone edges paired through a translation and optional traversal reversal.
     */
    struct CPeriodicMarker
    {
      std::size_t mIndexZoneI;  ///< Index of the first paired zone.
      std::size_t mIndexZoneJ;  ///< Index of the second paired zone.

      EFaceLocation mFaceLocationI;  ///< Local edge in the first paired zone.
      EFaceLocation mFaceLocationJ;  ///< Local edge in the second paired zone.

      std::string mNameI;  ///< Name of the first paired marker.
      std::string mNameJ;  ///< Name of the second paired marker.
      
      bool mIsReversed; ///< Whether neighbor traversal is reversed relative to owner traversal.

      std::array<as3double, 2> mTranslationVector; ///< Translation from the owner/first side to the neighbor/second side.
    

      /*!
       * @brief Getter function which returns the logical face direction of the first paired edge.
       *
       * @return logical face direction of the first paired edge.
       */
      ETypeFace GetiTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationI);
      }

      /*!
       * @brief Getter function which returns the logical face direction of the second paired edge.
       *
       * @return logical face direction of the second paired edge.
       */
      ETypeFace GetjTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationJ);
      }
    };

    /*!
     * @brief Constructor that initializes CPeriodicFamilyMarker from the specified data.
     *
     * @param[in] multizone_geometry_container geometry of all grid zones.
     * @param[in] iexternal_marker external marker on the first periodic side.
     * @param[in] jexternal_marker external marker on the second periodic side.
     * @param[in] translation translation from the first periodic side to the second.
     */
    CPeriodicFamilyMarker(const CMultizoneGeometry       *multizone_geometry_container,
                          const CExternalFamilyMarker    &iexternal_marker,
                          const CExternalFamilyMarker    &jexternal_marker,
                          const std::array<as3double, 2> &translation);

    /*!
     * @brief Getter function which returns the number of marker records in the family.
     *
     * @return number of marker records in the family.
     */
    std::size_t GetnMarkers(void) const {return mMarkers.size();}

    /*!
     * @brief Getter function which returns the marker records in the family.
     *
     * @return marker records in the family.
     */
    const auto& GetMarkers(void) const {return mMarkers;}

    /*!
     * @brief Getter function which returns the marker records in the family.
     *
     * @return marker records in the family.
     */
    auto& GetMarkers(void) {return mMarkers;}

    /*!
     * @brief Getter function which returns one family-local marker record.
     *
     * @param[in] iMarker zero-based marker index within its family.
     *
     * @return one family-local marker record.
     */
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}
  
  private:
    
    as3vector1d<CPeriodicMarker> mMarkers; ///< Zone-edge marker records in this family.
};



/*!
 * @brief Pairs of zone-edge markers connected by grid interfaces.
 */
class CInternalFamilyMarker
{
  public:
    

    /*!
     * @brief Two conforming zone edges paired with an optional traversal reversal.
     */
    struct CInternalMarker
    {
      std::size_t mIndexZoneI;  ///< Index of the first paired zone.
      std::size_t mIndexZoneJ;  ///< Index of the second paired zone.

      EFaceLocation mFaceLocationI;  ///< Local edge in the first paired zone.
      EFaceLocation mFaceLocationJ;  ///< Local edge in the second paired zone.

      bool mIsReversed; ///< Whether neighbor traversal is reversed relative to owner traversal.

      /*!
       * @brief Getter function which returns the logical face direction of the first paired edge.
       *
       * @return logical face direction of the first paired edge.
       */
      ETypeFace GetiTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationI);
      }

      /*!
       * @brief Getter function which returns the logical face direction of the second paired edge.
       *
       * @return logical face direction of the second paired edge.
       */
      ETypeFace GetjTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocationJ);
      }
    };

    /*!
     * @brief Constructor that initializes CInternalFamilyMarker from the specified data.
     *
     * @param[in] nMarker capacity reserved for marker records.
     */
    CInternalFamilyMarker(std::size_t nMarker)
    {
      mMarkers.reserve( nMarker );
    }

    /*!
     * @brief Getter function which returns the number of marker records in the family.
     *
     * @return number of marker records in the family.
     */
    std::size_t GetnMarkers(void) const {return mMarkers.size();}

    /*!
     * @brief Getter function which returns the marker records in the family.
     *
     * @return marker records in the family.
     */
    const auto& GetMarkers(void) const {return mMarkers;}

    /*!
     * @brief Getter function which returns the marker records in the family.
     *
     * @return marker records in the family.
     */
    auto& GetMarkers(void) {return mMarkers;}

    /*!
     * @brief Getter function which returns one family-local marker record.
     *
     * @param[in] iMarker zero-based marker index within its family.
     *
     * @return one family-local marker record.
     */
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}
  
  private:
    
    as3vector1d<CInternalMarker> mMarkers; ///< Zone-edge marker records in this family.
};


/*!
 * @brief Named physical boundary markers on zone edges.
 */
class CExternalFamilyMarker
{
  public:

    /*!
     * @brief One zone edge assigned to a named physical boundary.
     */
    struct CExternalMarker
    {
      std::size_t   mIndexZone;     ///< Index of the containing zone.
      EFaceLocation mFaceLocation;  ///< Local zone edge: IMIN, IMAX, JMIN or JMAX.

      /*!
       * @brief Getter function which returns the logical face direction associated with the local edge.
       *
       * @return logical face direction associated with the local edge.
       */
      ETypeFace GetTypeFace(void) const
      {
        return NMarkerUtility::GetTypeFace(mFaceLocation);
      }
    };

    /*!
     * @brief Constructor that initializes CExternalFamilyMarker from the specified data.
     *
     * @param[in] name name identifying the boundary marker.
     * @param[in] nMarker capacity reserved for marker records.
     */
    CExternalFamilyMarker(std::string name,
                          std::size_t nMarker)
      : mName(std::move(name))
    {
      mMarkers.reserve( nMarker );
    }

    /*!
     * @brief Getter function which returns the number of marker records in the family.
     *
     * @return number of marker records in the family.
     */
    std::size_t GetnMarkers(void) const {return mMarkers.size();}

    /*!
     * @brief Getter function which returns the boundary marker name.
     *
     * @return boundary marker name.
     */
    const std::string& GetName(void) const {return mName;}

    /*!
     * @brief Getter function which returns the marker records in the family.
     *
     * @return marker records in the family.
     */
    const auto& GetMarkers(void) const {return mMarkers;}

    /*!
     * @brief Getter function which returns the marker records in the family.
     *
     * @return marker records in the family.
     */
    auto& GetMarkers(void) {return mMarkers;}

    /*!
     * @brief Getter function which returns one family-local marker record.
     *
     * @param[in] iMarker zero-based marker index within its family.
     *
     * @return one family-local marker record.
     */
    const auto& GetMarker(std::size_t iMarker) const {return mMarkers[iMarker];}

  private:
    
    std::string mName;                      ///< Name identifying this boundary marker.
    as3vector1d<CExternalMarker> mMarkers;  ///< Zone-edge marker records in this family.
};




