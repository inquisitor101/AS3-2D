#pragma once

/*!
 * @file
 * @brief Face connectivity grouped by owner-side logical direction.
 */


#include <cstddef>

#include "geometry/face/CBoundaryFaceBatch.hpp"
#include "geometry/face/CInterfaceFaceBatch.hpp"
#include "geometry/face/CInteriorFaceBatch.hpp"
#include "geometry/face/CFaceCollection.hpp"
#include "geometry/face/CFlattenedFaceIndex.hpp"

class CMultizoneGeometry;


/*!
 * @brief Interior, boundary and interface face collections in one owner-side logical direction.
 */
class CDirectionalFaceConnectivity
{

    /*!
     * @brief NOTE, this class stores faces across all zones in one logical direction:
     *   [interior][boundary][interface]. An interface uses its owner's direction.
     *
     * @param[in] face_direction logical I- or J-face direction.
     */
	public:

    explicit CDirectionalFaceConnectivity(ETypeFace face_direction) : mFaceDirection(face_direction) {}

		/*!
		 * @brief Function that constructs face records from zone-edge connectivity.
		 *
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		void InitializeFaces(const CMultizoneGeometry *multizone_geometry_container);

    /*!
     * @brief Function that resolves a direction-local face index into its category, batch and batch-local face
     *   index.
     *
     * @param[in] i zero-based flattened index local to this collection or direction.
     *
     * @return face category, category-local batch index and batch-local face index.
     */
    CFlattenedFaceIndex ResolveFace(std::size_t i) const;

    /*!
     * @brief Getter function which returns the interior-face collection.
     *
     * @return interior-face collection.
     */
    const auto& GetInteriorFaces(void)  const {return mInteriorFaces;}

    /*!
     * @brief Getter function which returns the physical boundary-face collection.
     *
     * @return physical boundary-face collection.
     */
    const auto& GetBoundaryFaces(void)  const {return mBoundaryFaces;}

    /*!
     * @brief Getter function which returns the interface-face collection.
     *
     * @return interface-face collection.
     */
    const auto& GetInterfaceFaces(void) const {return mInterfaceFaces;}

		/*!
		 * @brief Getter function which returns the number of interior faces.
		 *
		 * @return number of interior faces.
		 */
		std::size_t GetnInteriorFaces(void)  const { return mInteriorFaces.GetnFace();  }

		/*!
		 * @brief Getter function which returns the number of physical boundary faces.
		 *
		 * @return number of physical boundary faces.
		 */
		std::size_t GetnBoundaryFaces(void)  const { return mBoundaryFaces.GetnFace();  }

		/*!
		 * @brief Getter function which returns the number of interface faces.
		 *
		 * @return number of interface faces.
		 */
		std::size_t GetnInterfaceFaces(void) const { return mInterfaceFaces.GetnFace(); }

		/*!
		 * @brief Getter function which returns the total number of faces across the three categories.
		 *
		 * @return total number of faces across the three categories.
		 */
		std::size_t GetnFacesTotal(void) const
		{
			return GetnInteriorFaces() + GetnBoundaryFaces() + GetnInterfaceFaces();
		}

    /*!
     * @brief Getter function which returns the logical I- or J-face direction.
     *
     * @return logical I- or J-face direction.
     */
    ETypeFace GetFaceDirection(void) const { return mFaceDirection; }

	private:
    ETypeFace mFaceDirection; ///< Logical I- or J-face direction.

    CFaceCollection<CInteriorFaceBatch>  mInteriorFaces;   ///< Interior faces in the selected logical direction.
    CFaceCollection<CBoundaryFaceBatch>  mBoundaryFaces;   ///< Physical boundary faces in the selected logical direction.
    CFaceCollection<CInterfaceFaceBatch> mInterfaceFaces;  ///< Interface faces grouped by their owner-side logical direction.


		/*!
		 * @brief Function that collects same-zone interior faces in this logical direction.
		 *
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		void InitializeInteriorFaces(const CMultizoneGeometry *multizone_geometry_container);

		/*!
		 * @brief Function that collects physical boundary faces in this logical direction.
		 *
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		void InitializeBoundaryFaces(const CMultizoneGeometry *multizone_geometry_container);

		/*!
		 * @brief Function that collects interfaces according to their owner-side logical direction.
		 *
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		void InitializeInterfaceFaces(const CMultizoneGeometry *multizone_geometry_container);
};
