#pragma once

/*!
 * @file
 * @brief Paired face records and batches between two zone edges.
 */


#include <array>
#include <cstddef>
#include <string>

#include "common/AS3Options.hpp"
#include "MarkerFamilies.hpp"

class CMultizoneGeometry;


/*!
 * @brief A pair of zone-local elements connected by an interface face.
 */
struct CInterfaceFace
{
  std::size_t mIndexElementOwner;     ///< Element index local to the owner zone.
  std::size_t mIndexElementNeighbor;  ///< Element index local to the neighbor zone.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Paired zone-edge faces with owner, neighbor and orientation metadata.
 */
class CInterfaceFaceBatch
{

  /*!
   * @brief Vector storing batch-local face records.
   */
  using AInterfaceFaceVector = as3vector1d<CInterfaceFace>;
  
  public:

    /*!
     * @brief Constructor that initializes CInterfaceFaceBatch from the specified data.
     *
     * @param[in] multizone_geometry_container geometry of all grid zones.
     * @param[in] marker zone-edge marker defining this batch.
     */
    CInterfaceFaceBatch(const CMultizoneGeometry                     *multizone_geometry_container,
                        const CInternalFamilyMarker::CInternalMarker  &marker);

    /*!
     * @brief Constructor that initializes CInterfaceFaceBatch from the specified data.
     *
     * @param[in] multizone_geometry_container geometry of all grid zones.
     * @param[in] marker zone-edge marker defining this batch.
     */
    CInterfaceFaceBatch(const CMultizoneGeometry                     *multizone_geometry_container,
                        const CPeriodicFamilyMarker::CPeriodicMarker  &marker);

    /*!
     * @brief Getter function which returns the face records in this batch.
     *
     * @return face records in this batch.
     */
    const AInterfaceFaceVector& GetFaces(void) const
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
    const CInterfaceFace& GetFace(std::size_t i) const
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
     * @brief Getter function which returns the owner-zone index.
     *
     * @return owner-zone index.
     */
    std::size_t GetOwnerZoneIndex(void) const {return mIndexZoneOwner;}

    /*!
     * @brief Getter function which returns the neighbor-zone index.
     *
     * @return neighbor-zone index.
     */
    std::size_t GetNeighborZoneIndex(void) const {return mIndexZoneNeighbor;}

    /*!
     * @brief Getter function which returns the whether the neighbor traversal is reversed relative to the
     *   owner.
     *
     * @return whether the neighbor traversal is reversed relative to the owner.
     */
    bool IsNeighborReversed(void) const {return mIsReversed;}

    /*!
     * @brief Getter function which returns the owner-side local edge.
     *
     * @return owner-side local edge.
     */
    EFaceLocation GetOwnerFaceLocation(void) const {return mFaceLocationOwner;}

    /*!
     * @brief Getter function which returns the neighbor-side local edge.
     *
     * @return neighbor-side local edge.
     */
    EFaceLocation GetNeighborFaceLocation(void) const {return mFaceLocationNeighbor;}

		/*!
		 * @brief Getter function which returns the whether the owner edge is IMAX or JMAX.
		 *
		 * @return whether the owner edge is IMAX or JMAX.
		 */
		bool IsOwnerMaxFace(void) const {return mIsMaxFaceOwner;}

		/*!
		 * @brief Getter function which returns the whether the neighbor edge is IMAX or JMAX.
		 *
		 * @return whether the neighbor edge is IMAX or JMAX.
		 */
		bool IsNeighborMaxFace(void) const {return mIsMaxFaceNeighbor;}

    /*!
     * @brief Getter function which returns the owner-side marker name.
     *
     * @return owner-side marker name.
     */
    const std::string& GetOwnerFaceName() const { return mFaceNameOwner; }

    /*!
     * @brief Getter function which returns the neighbor-side marker name.
     *
     * @return neighbor-side marker name.
     */
    const std::string& GetNeighborFaceName() const { return mFaceNameNeighbor; }

    /*!
     * @brief Getter function which returns the owner-side logical face direction.
     *
     * @return owner-side logical face direction.
     */
    ETypeFace GetOwnerFaceDirection(void) const { return mFaceDirectionOwner; }

    /*!
     * @brief Getter function which returns the neighbor-side logical face direction.
     *
     * @return neighbor-side logical face direction.
     */
    ETypeFace GetNeighborFaceDirection(void) const { return mFaceDirectionNeighbor; }

  private:
    std::size_t mIndexZoneOwner;     ///< Owner-zone index.
    std::size_t mIndexZoneNeighbor;  ///< Neighbor-zone index.

    ETypeFace mFaceDirectionOwner;     ///< Logical face direction in the owner zone.
    ETypeFace mFaceDirectionNeighbor;  ///< Logical face direction in the neighbor zone.

    EFaceLocation mFaceLocationOwner;     ///< Local edge in the owner zone.
    EFaceLocation mFaceLocationNeighbor;  ///< Local edge in the neighbor zone.

    std::string mFaceNameOwner;     ///< Name of the owner-side marker.
    std::string mFaceNameNeighbor;  ///< Name of the neighbor-side marker.

		bool mIsMaxFaceOwner;     ///< Whether the owner edge is IMAX or JMAX.
		bool mIsMaxFaceNeighbor;  ///< Whether the neighbor edge is IMAX or JMAX.

    bool mIsReversed; ///< Whether neighbor traversal is reversed relative to owner traversal.

    std::array<as3double, 2> mTranslationVector; ///< Translation from the owner/first side to the neighbor/second side.

    AInterfaceFaceVector mFaces; ///< Face records local to this batch.

    /*!
     * @brief Function that constructs face records from zone-edge connectivity.
     *
     * @param[in] multizone_geometry_container geometry of all grid zones.
     */
    void InitializeFaces(const CMultizoneGeometry *multizone_geometry_container);

		/*!
		 * @brief Function that checks paired face coordinates, orientation and periodic translation.
		 *
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		void CheckConformity(const CMultizoneGeometry *multizone_geometry_container) const;
};
