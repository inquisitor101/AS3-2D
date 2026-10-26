#pragma once 

/*!
 * @file
 * @brief Interior, boundary and interface face connectivity and batch-local indexing.
 */


#include <array>
#include <cstddef>
#include <string>
#include <utility>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "marker_structure.hpp"
#include "flattened_index_structure.hpp"


// Forward declaration to avoid compiler problems.
class CMultizoneGeometry;
class CSinglezoneGeometry;




/*!
 * @brief The batch index and batch-local face index within one face collection.
 */
struct CFaceAddress
{
  std::size_t mIndexBatch;  ///< Batch index local to the selected face category.
  std::size_t mIndexFace;   ///< Face index local to its batch.
};




//-----------------------------------------------------------------------------------


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

//-----------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------

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

//-----------------------------------------------------------------------------------


/*!
 * @brief A collection of batches with a collection-local flattened face index.
 *
 * @tparam TBatch face batch type stored by the collection.
 */
template<typename TBatch>
class CFaceCollection
{
  private:

    /*!
     * @brief A half-open range of flattened face indices assigned to one batch.
     */
    struct CFaceRange
    {

      /*!
       * @brief Function that tests membership in this half-open flattened face range.
       *
       * @param[in] i flattened index tested against this half-open range.
       *
       * @return true if the index lies in [mBegin,mEnd); false otherwise.
       */
      bool Contains(std::size_t i) const { return i >= mBegin && i < mEnd; }

      std::size_t mBegin;  ///< First flattened face index in the range.
      std::size_t mEnd;    ///< One-past-last flattened face index in the range.
    };

  public:

    /*!
     * @brief Default constructor that initializes an empty CFaceCollection.
     */
    CFaceCollection(void) = default;

    /*!
     * @brief Function that takes ownership of face batches and rebuilds their flattened ranges.
     *
     * @param[in,out] batches face batches whose ownership is transferred to this collection.
     *
     * @note Initialization replaces the batches and their ranges; complete it before parallel traversal.
     */
    void Initialize(as3vector1d<TBatch>&& batches)
    {
      mBatches = std::move(batches);

      mBatchRanges.clear();
      mBatchRanges.reserve( mBatches.size() );

      mNFace = 0;
      for(const auto& batch : mBatches)
      {
        const std::size_t nFace = batch.GetnFace();
        mBatchRanges.emplace_back( CFaceRange{mNFace, mNFace+nFace} );
        mNFace += nFace;
      }
    }

    /*!
     * @brief Function that resolves a collection-local flattened face index into a batch and batch-local index.
     *
     * @param[in] i zero-based flattened index local to this collection or direction.
     *
     * @return batch index and batch-local face index; terminates if the index is invalid.
     */
    CFaceAddress ResolveFace(std::size_t i) const
    {
      // The flat index is local to this face collection.
      for(std::size_t iBatch=0; iBatch<mBatchRanges.size(); iBatch++)
      {
        const CFaceRange& range = mBatchRanges[iBatch];
        if( range.Contains(i) ) return CFaceAddress{ iBatch, i - range.mBegin };
      }
      ERROR("Invalid collection-local face index.");
    }

    /*!
     * @brief Getter function which returns the number of faces.
     *
     * @return number of faces.
     */
    std::size_t GetnFace(void) const { return mNFace; }

    /*!
     * @brief Getter function which returns the number of face batches.
     *
     * @return number of face batches.
     */
    std::size_t GetnBatch(void) const { return mBatches.size(); }

    /*!
     * @brief Getter function which returns the face batches in this collection.
     *
     * @return face batches in this collection.
     */
    const as3vector1d<TBatch>& GetBatches(void) const
    {
      return mBatches;
    }

    /*!
     * @brief Getter function which returns one batch selected by a collection-local index.
     *
     * @param[in] iBatch batch index local to the selected face category.
     *
     * @return one batch selected by a collection-local index.
     */
    const TBatch& GetBatch(std::size_t iBatch) const
    {
#if DEBUG
      if( iBatch >= mBatches.size() ) ERROR("Invalid face batch index.");
#endif
      return mBatches[iBatch];
    }

  private:
    std::size_t mNFace = 0; ///< Number of faces in this collection or interface.

    as3vector1d<TBatch>     mBatches;      ///< Owned face batches.
    as3vector1d<CFaceRange> mBatchRanges;  ///< Collection-local flattened ranges corresponding to the batches.
};

//-----------------------------------------------------------------------------------


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





