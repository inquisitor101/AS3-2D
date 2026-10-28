#pragma once 

/*!
 * @file
 * @brief Multizone topology, zone grids and element coordinate storage.
 */


#include <cstddef>
#include <exception>
#include <memory>
#include <string>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "openmp_structure.hpp"
#include "input_structure.hpp"
#include "marker_structure.hpp"
#include "face_structure.hpp"
#include "load_balance_faces.hpp"
#include "flattened_index_structure.hpp"

// Forward declaration to avoid compiler problems.
class CSinglezoneGeometry;
class CElementGeometry;


/*!
 * @brief Geometry and flattened element/face traversal for the complete multizone grid.
 */
class CMultizoneGeometry
{
	public:

		/*!
		 * @brief Constructor of CMultizoneGeometry, which is responsible for the entire grid geometry.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 */
		CMultizoneGeometry(const CConfig *config_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CMultizoneGeometry(void);


		/*!
		 * @brief Getter function which returns a vector of grid zone containers.
		 *
		 * @return mSinglezoneGeometry
		 */
		const as3vector1d<std::unique_ptr<CSinglezoneGeometry>>& GetSinglezoneGeometry(void) const
    {
      return mSinglezoneGeometry;
    }

		/*!
		 * @brief Getter function which returns a specific grid zone container.
		 *
		 * @param[in] iZone zero-based zone index.
		 *
		 * @return mSinglezoneGeometry[iZone]
		 */
		CSinglezoneGeometry* GetSinglezoneGeometry(std::size_t iZone) const { return mSinglezoneGeometry[iZone].get(); }

    /*!
     * @brief Function that builds flattened element and face indices and face load-balancing permutations.
     *
     * @param[in] openmp_container OpenMP team configuration.
     */
    void InitializeGridTopology(const COpenMP *openmp_container);

		/*!
		 * @brief Getter function which returns the connectivity classified by the owner-side I
		 *   direction.
		 *
		 * @return connectivity classified by the owner-side I direction.
		 */
		const CDirectionalFaceConnectivity& GetMultizoneIFaces(void) const { return mMultizoneIFaces; }

		/*!
		 * @brief Getter function which returns the connectivity classified by the owner-side J
		 *   direction.
		 *
		 * @return connectivity classified by the owner-side J direction.
		 */
		const CDirectionalFaceConnectivity& GetMultizoneJFaces(void) const { return mMultizoneJFaces; }

    /*!
     * @brief Getter function which returns the number of zones.
     *
     * @return number of zones.
     */
    std::size_t GetnZone(void) const { return mNZone; }

		/*!
		 * @brief Getter function which returns the total number of owner-classified I-faces.
		 *
		 * @return total number of owner-classified I-faces.
		 */
		std::size_t GetnIFace(void) const { return mMultizoneIFaces.GetnFacesTotal(); }

		/*!
		 * @brief Getter function which returns the total number of owner-classified J-faces.
		 *
		 * @return total number of owner-classified J-faces.
		 */
		std::size_t GetnJFace(void) const { return mMultizoneJFaces.GetnFacesTotal(); }

    /*!
     * @brief Getter function which returns the total number of elements across all zones.
     *
     * @return total number of elements across all zones.
     */
    std::size_t GetnElemTotal(void) const { return mNElemTotal; }

		/*!
		 * @brief Getter function which returns the zone-local element address for a global traversal
		 *   index.
		 *
		 * @param[in] i zero-based traversal index.
		 *
		 * @return zone-local element address for a global traversal index.
		 */
		CFlattenedElementIndex GetFlattenedIndexVolumeElement(std::size_t i) const
		{
#if DEBUG
			if( i >= mNElemTotal ) ERROR("Index is out of range.");
#endif
			return mFlattenedIndexVolumeElement[i];
		}

    /*!
     * @brief Getter function which returns the category-local I-face address before load balancing.
     *
     * @param[in] i zero-based traversal index.
     *
     * @return category-local I-face address before load balancing.
     */
    CFlattenedFaceIndex GetFlattenedIndexIFace(std::size_t i) const
    {
#if DEBUG
      if( i >= GetnIFace() ) ERROR("Index is out of range.");
#endif
      return mFlattenedIndexIFace[i];
    }

    /*!
     * @brief Getter function which returns the category-local J-face address before load balancing.
     *
     * @param[in] i zero-based traversal index.
     *
     * @return category-local J-face address before load balancing.
     */
    CFlattenedFaceIndex GetFlattenedIndexJFace(std::size_t i) const
    {
#if DEBUG
      if( i >= GetnJFace() ) ERROR("Index is out of range.");
#endif
      return mFlattenedIndexJFace[i];
    }

    /*!
     * @brief Getter function which returns the I-face load-balancing permutation.
     *
     * @return I-face load-balancing permutation.
     */
    const CLoadBalancedFacePermutation& GetIFaceLoadBalancedPermutation(void) const
    {
      return mIFaceLoadBalancedPermutation;
    }

    /*!
     * @brief Getter function which returns the J-face load-balancing permutation.
     *
     * @return J-face load-balancing permutation.
     */
    const CLoadBalancedFacePermutation& GetJFaceLoadBalancedPermutation(void) const
    {
      return mJFaceLoadBalancedPermutation;
    }

    /*!
     * @brief Getter function which returns the category-local face address at the permuted I-face position.
     *
     * @param[in] i zero-based traversal index.
     *
     * @return category-local face address at the permuted I-face position.
     */
    CFlattenedFaceIndex GetFlattenedIndexIFaceLoadBalanced(std::size_t i) const
    {
#if DEBUG
      if( i >= mIFaceLoadBalancedPermutation.GetnFace() ) ERROR("Index is out of range.");
#endif
      return mFlattenedIndexIFace[ mIFaceLoadBalancedPermutation.GetPermutationIndex(i) ];
    }

    /*!
     * @brief Getter function which returns the category-local face address at the permuted J-face position.
     *
     * @param[in] i zero-based traversal index.
     *
     * @return category-local face address at the permuted J-face position.
     */
    CFlattenedFaceIndex GetFlattenedIndexJFaceLoadBalanced(std::size_t i) const
    {
#if DEBUG
      if( i >= mJFaceLoadBalancedPermutation.GetnFace() ) ERROR("Index is out of range.");
#endif
      return mFlattenedIndexJFace[ mJFaceLoadBalancedPermutation.GetPermutationIndex(i) ];
    }

    /*!
     * @brief Getter function which returns the interface batches classified by the owner-side I direction.
     *
     * @return interface batches classified by the owner-side I direction.
     */
    const auto& GetInterfaceBatchesIFace(void) const
    {
      return mMultizoneIFaces.GetInterfaceFaces().GetBatches();
    }

    /*!
     * @brief Getter function which returns the interface batches classified by the owner-side J direction.
     *
     * @return interface batches classified by the owner-side J direction.
     */
    const auto& GetInterfaceBatchesJFace(void) const
    {
      return mMultizoneJFaces.GetInterfaceFaces().GetBatches();
    } 

    /*!
     * @brief Getter function which returns the physical boundary batches in the I direction.
     *
     * @return physical boundary batches in the I direction.
     */
    const auto& GetBoundaryBatchesIFace(void) const
    {
      return mMultizoneIFaces.GetBoundaryFaces().GetBatches();
    }

    /*!
     * @brief Getter function which returns the physical boundary batches in the J direction.
     *
     * @return physical boundary batches in the J direction.
     */
    const auto& GetBoundaryBatchesJFace(void) const
    {
      return mMultizoneJFaces.GetBoundaryFaces().GetBatches();
    }

    /*!
     * @brief Getter function which returns the internal marker families.
     *
     * @return internal marker families.
     */
    const auto& GetInternalFamilyMarkers(void) const { return mInternalFamilyMarkers; }

    /*!
     * @brief Getter function which returns the physical boundary marker families.
     *
     * @return physical boundary marker families.
     */
    const auto& GetExternalFamilyMarkers(void) const { return mExternalFamilyMarkers; }

    /*!
     * @brief Getter function which returns the periodic marker families.
     *
     * @return periodic marker families.
     */
    const auto& GetPeriodicFamilyMarkers(void) const { return mPeriodicFamilyMarkers; }

	protected:

	private:
		std::size_t mNZone = 0;                                                 ///< Number of grid zones.
		std::size_t mNElemTotal = 0;                                            ///< Total number of elements across all zones.
		as3vector1d<std::unique_ptr<CSinglezoneGeometry>> mSinglezoneGeometry;  ///< Owned geometry containers indexed by zone.

    as3vector1d<CInternalFamilyMarker> mInternalFamilyMarkers;  ///< Internal zone-edge marker families.
    as3vector1d<CExternalFamilyMarker> mExternalFamilyMarkers;  ///< Named physical boundary marker families.
    as3vector1d<CPeriodicFamilyMarker> mPeriodicFamilyMarkers;  ///< Periodic zone-edge marker families.

		as3vector1d<CFlattenedElementIndex> mFlattenedIndexVolumeElement;  ///< Global element traversal mapped to zone-local element addresses.
    as3vector1d<CFlattenedFaceIndex>    mFlattenedIndexIFace;                      ///< I-direction traversal mapped to category-local batch and face addresses.
    as3vector1d<CFlattenedFaceIndex>    mFlattenedIndexJFace;                      ///< J-direction traversal mapped to category-local batch and face addresses.

    CLoadBalancedFacePermutation mIFaceLoadBalancedPermutation;  ///< Cost-weighted permutation of I-direction face indices.
    CLoadBalancedFacePermutation mJFaceLoadBalancedPermutation;  ///< Cost-weighted permutation of J-direction face indices.

    CDirectionalFaceConnectivity mMultizoneIFaces{ETypeFace::IFACE};  ///< Face connectivity classified by the owner-side I direction.
    CDirectionalFaceConnectivity mMultizoneJFaces{ETypeFace::JFACE};  ///< Face connectivity classified by the owner-side J direction.

		

		/*!
		 * @brief Function that checks the existence of the grid files and reports the output.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 */
		void CheckExistanceGridFiles(const CConfig *config_container);


    /*!
     * @brief Function that imports AS3 grid metadata, zone coordinates and zone-edge markers.
     *
     * @param[in] config_container configuration/dictionary container.
     */
    void ImportGrid(const CConfig *config_container);


		// Disable default constructor.
		CMultizoneGeometry(void) = delete;
		// Disable default copy constructor.
		CMultizoneGeometry(const CMultizoneGeometry&) = delete;
		// Disable default copy operator.
		CMultizoneGeometry& operator=(CMultizoneGeometry&) = delete;	
};

//-----------------------------------------------------------------------------------

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
    CSinglezoneGeometry(const std::string &gridfile, 
                        std::size_t        iZone, 
                        std::size_t        nPoly, 
                        std::size_t        niElem, 
                        std::size_t        njElem, 
                        ETypeDOF           nodalDistribution, 
                        bool               isAffine);


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
    ETypeDOF GetTypeNodalDistribution(void) const {return mNodalDistribution;}

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

    ETypeDOF    mNodalDistribution;                                   ///< Distribution of geometry nodes within each element.
		as3vector1d<std::unique_ptr<CElementGeometry>> mElementGeometry;  ///< Container with the element geometry.


		// Disable default constructor.
		CSinglezoneGeometry(void) = delete;
		// Disable default copy constructor.
		CSinglezoneGeometry(const CSinglezoneGeometry&) = delete;
		// Disable default copy operator.
		CSinglezoneGeometry& operator=(CSinglezoneGeometry&) = delete;	
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Coordinates of the geometry nodes of one quadrilateral element.
 */
class CElementGeometry
{
	public:

    /*!
     * @brief Constructor that initializes CElementGeometry from the specified data.
     *
     * @param[in] nNode2D number of geometry nodes in the element.
     */
    explicit CElementGeometry(std::size_t nNode2D);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CElementGeometry(void);
	

		/*!
		 * @brief Getter function which returns the coordinates at the solution DOFs.
		 *
		 * @return mCoordSolDOFs.
		 */
		CMatrixAS3<as3double> &GetCoordSolDOFs(void) {return mCoordSolDOFs;}

	protected:

	private:
		CMatrixAS3<as3double> mCoordSolDOFs;  ///< Coordinates at the solution DOFs.


		// Disable default constructor.
		CElementGeometry(void) = delete;
		// Disable default copy constructor.
		CElementGeometry(const CElementGeometry&) = delete;
		// Disable default copy operator.
		CElementGeometry& operator=(CElementGeometry&) = delete;	
};







