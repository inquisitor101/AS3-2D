#pragma once

/*!
 * @file
 * @brief Physical boundary face records and batches for one zone edge.
 */


#include <cstddef>
#include <string>

#include "common/AS3Options.hpp"
#include "MarkerFamilies.hpp"

class CSinglezoneGeometry;


/*!
 * @brief The adjacent zone-local element of a physical boundary face.
 */
struct CBoundaryFace
{
  std::size_t mIndexElement; ///< Zone-local element adjacent to this boundary face.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Physical boundary faces sharing one zone edge and marker name.
 */
class CBoundaryFaceBatch
{

  /*!
   * @brief Vector storing batch-local face records.
   */
  using ABoundaryFaceVector = as3vector1d<CBoundaryFace>;
  
  public:

    /*!
     * @brief Constructor that initializes CBoundaryFaceBatch from the specified data.
     *
     * @param[in] singlezone_geometry_container geometry of the current zone.
     * @param[in] marker zone-edge marker defining this batch.
     * @param[in] name name identifying the boundary marker.
     */
    CBoundaryFaceBatch(const CSinglezoneGeometry                    *singlezone_geometry_container,
                       const CExternalFamilyMarker::CExternalMarker &marker,
                       const std::string                            &name);

    /*!
     * @brief Getter function which returns the face records in this batch.
     *
     * @return face records in this batch.
     */
    const ABoundaryFaceVector& GetFaces(void) const
    {
      return mFaces;
    }

    /*!
     * @brief Getter function which returns one batch-local face record.
     *
     * @param[in] i zero-based face index local to this batch.
     *
     * @return one batch-local face record.
     */
    const CBoundaryFace& GetFace(std::size_t i) const
    {
      return mFaces[i];
    }

    /*!
     * @brief Getter function which returns the number of faces.
     *
     * @return number of faces.
     */
    std::size_t GetnFace(void) const {return mFaces.size();}

    /*!
     * @brief Getter function which returns the zone index owning this batch.
     *
     * @return zone index owning this batch.
     */
    std::size_t GetZoneIndex(void) const {return mIndexZone;}

    /*!
     * @brief Getter function which returns the logical I- or J-face direction.
     *
     * @return logical I- or J-face direction.
     */
    ETypeFace GetFaceDirection(void) const {return mFaceDirection;}

    /*!
     * @brief Getter function which returns the boundary marker name.
     *
     * @return boundary marker name.
     */
    std::string GetName(void) const {return mName;}

    /*!
     * @brief Getter function which returns the local zone-edge location.
     *
     * @return local zone-edge location.
     */
    EFaceLocation GetFaceLocation(void) const {return mFaceLocation;}

    /*!
     * @brief Getter function which returns the whether this edge is IMAX or JMAX.
     *
     * @return whether this edge is IMAX or JMAX.
     */
    bool IsMaxFace(void) const {return mIsMaxFace;}

  private:
    ETypeFace           mFaceDirection;  ///< Logical I- or J-face direction.
    EFaceLocation       mFaceLocation;   ///< Local zone edge: IMIN, IMAX, JMIN or JMAX.
    std::size_t         mIndexZone;      ///< Index of the containing zone.
    std::string         mName;           ///< Name identifying this boundary marker.
    bool                mIsMaxFace;      ///< Whether this edge is IMAX or JMAX.
    ABoundaryFaceVector mFaces;          ///< Face records local to this batch.

    /*!
     * @brief Function that constructs face records from zone-edge connectivity.
     *
     * @param[in] singlezone_geometry_container geometry of the current zone.
     */
    void InitializeFaces(const CSinglezoneGeometry *singlezone_geometry_container);
};
