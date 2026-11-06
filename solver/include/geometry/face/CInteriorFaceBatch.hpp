#pragma once

/*!
 * @file
 * @brief Interior face records and batches for one zone and logical direction.
 */


#include <cstddef>

#include "common/AS3Options.hpp"

class CSinglezoneGeometry;


/*!
 * @brief The minus-side element of an interior face within a zone.
 */
struct CInteriorFace
{
  std::size_t mIndexElementMinus; ///< Zone-local minus-side element; add the batch offset to obtain the plus-side element.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Interior faces sharing one zone, logical direction and neighbor offset.
 */
class CInteriorFaceBatch
{

  /*!
   * @brief Vector storing batch-local face records.
   */
  using AInteriorFaceVector = as3vector1d<CInteriorFace>;
  
  public:

    /*!
     * @brief Constructor that initializes CInteriorFaceBatch from the specified data.
     *
     * @param[in] face_direction logical I- or J-face direction.
     * @param[in] singlezone_geometry_container geometry of the current zone.
     */
    CInteriorFaceBatch(ETypeFace                  face_direction, 
                       const CSinglezoneGeometry *singlezone_geometry_container);

    /*!
     * @brief Getter function which returns the face records in this batch.
     *
     * @return face records in this batch.
     */
    const AInteriorFaceVector& GetFaces(void) const
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
    const CInteriorFace& GetFace(std::size_t i) const
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
     * @brief Getter function which returns the offset from the minus-side element to its plus-side neighbor.
     *
     * @return offset from the minus-side element to its plus-side neighbor.
     */
    std::size_t GetNeighborElementOffset(void) const {return mOffsetElementPlus;}

  private:
    ETypeFace           mFaceDirection;      ///< Logical I- or J-face direction.
    std::size_t         mIndexZone;          ///< Index of the containing zone.
    std::size_t         mOffsetElementPlus;  ///< Offset from the minus-side element to its plus-side neighbor.
    AInteriorFaceVector mFaces;              ///< Face records local to this batch.

    /*!
     * @brief Function that constructs face records from zone-edge connectivity.
     *
     * @param[in] singlezone_geometry_container geometry of the current zone.
     */
    void InitializeFaces(const CSinglezoneGeometry *singlezone_geometry_container);
};
