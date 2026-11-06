#pragma once

/*!
 * @file
 * @brief Structured grid geometry for one zone.
 */


#include <cstddef>
#include <memory>
#include <string>

#include "common/AS3Options.hpp"

class CElementGeometry;


/*!
 * @brief A class used for storing a (single) zone grid.
 */
class CSinglezoneGeometry
{
	public:

    /*!
     * @brief Constructor that initializes CSinglezoneGeometry from the specified data.
     *
     * @param[in] gridfile associated grid filename.
     * @param[in] iZone zero-based zone index.
     * @param[in] nPoly polynomial order in each reference direction.
     * @param[in] niElem number of elements in the logical i-direction.
     * @param[in] njElem number of elements in the logical j-direction.
     * @param[in] nodalDistribution distribution of the geometry nodes.
     * @param[in] isAffine whether the grid declares affine element mappings.
     */
    CSinglezoneGeometry(const std::string      &gridfile,
                        std::size_t            iZone,
                        std::size_t            nPoly,
                        std::size_t            niElem,
                        std::size_t            njElem,
                        ETypeNodalDistribution nodalDistribution,
                        bool                   isAffine);


		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CSinglezoneGeometry(void);

    /*!
     * @brief Getter function which returns the value of mZoneID.
     *
     * @return mZoneID.
     */
    std::size_t GetZoneID(void) const {return mZoneID;}

		/*!
		 * @brief Getter function which returns the value of mGridFile.
		 *
		 * @return mGridFile.
		 */
		const std::string &GetGridFile(void) const {return mGridFile;}

		/*!
		 * @brief Getter function which returns the element of the specified index.
		 *
		 * @param[in] index zero-based entry index.
		 *
		 * @return mElementGeometry[index].
		 */
		CElementGeometry *GetElementGeometry(std::size_t index) const {return mElementGeometry[index].get();}

    /*!
     * @brief Getter function which returns the total number of elements in this zone.
     *
     * @return mElementGeometry.size().
     */
    std::size_t GetnElem(void) const {return mElementGeometry.size();}

    /*!
     * @brief Getter function which returns the number of elements in the logical i-direction.
     *
     * @return mNiElem.
     */
    std::size_t GetniElem(void) const {return mNiElem;}

    /*!
     * @brief Getter function which returns the number of elements in the logical j-direction.
     *
     * @return mNjElem.
     */
    std::size_t GetnjElem(void) const {return mNjElem;}

    /*!
     * @brief Getter function which returns the number of grid nodes in 2D.
     *
     * @return (mNPolyGrid+1)*(mNPolyGrid+1).
     */
    std::size_t GetnNodeGrid2D(void) const {return (mNPolyGrid+1)*(mNPolyGrid+1);}

    /*!
     * @brief Getter function which returns the polynomial order of the grid element.
     *
     * @return mNPolyGrid.
     */
    std::size_t GetnPolyGrid(void) const {return mNPolyGrid;}

    /*!
     * @brief Getter function which returns the number of spatial dimensions, which is two.
     *
     * @return number of spatial dimensions, which is two.
     */
    std::size_t GetnDim(void) const {return 2;}

    /*!
     * @brief Getter function which returns the distribution of the geometry nodes.
     *
     * @return distribution of the geometry nodes.
     */
    ETypeNodalDistribution GetTypeNodalDistribution(void) const {return mNodalDistribution;}

    /*!
     * @brief Function that returns zone-local element indices in increasing traversal order along a zone edge.
     *
     * @param[in] face_location local edge: IMIN, IMAX, JMIN or JMAX.
     *
     * @return zone-local element indices ordered along the selected edge.
     */
    as3vector1d<std::size_t> ComputeSurfaceElementIndices(EFaceLocation face_location) const;

    /*!
     * @brief Function that returns reference-element node indices along a local edge.
     *
     * @param[in] face_location local edge: IMIN, IMAX, JMIN or JMAX.
     *
     * @return reference-element node indices ordered along the selected edge.
     */
    as3vector1d<std::size_t> ComputeSurfaceNodeIndices(EFaceLocation face_location) const;


	protected:

	private:
    std::size_t mZoneID;           ///< Current zone index.
		std::string	mGridFile;         ///< Associated grid file name of this zone.
    std::size_t mNPolyGrid;        ///< Polynomial order of the grid element.
		std::size_t mNiElem = 0;       ///< Number of elements in i-direction.
    std::size_t mNjElem = 0;       ///< Number of elements in j-direction.	
    bool        mIsAffine = false; ///< Whether the imported grid declares affine element mappings.

    ETypeNodalDistribution mNodalDistribution;                                   ///< Distribution of geometry nodes within each element.
		as3vector1d<std::unique_ptr<CElementGeometry>> mElementGeometry;  ///< Container with the element geometry.


		// Disable default constructor.
		CSinglezoneGeometry(void) = delete;
		// Disable default copy constructor.
		CSinglezoneGeometry(const CSinglezoneGeometry&) = delete;
		// Disable default copy operator.
		CSinglezoneGeometry& operator=(CSinglezoneGeometry&) = delete;	
};
