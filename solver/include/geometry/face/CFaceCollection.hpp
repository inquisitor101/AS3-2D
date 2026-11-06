#pragma once

/*!
 * @file
 * @brief Batch collections and their flattened face addresses.
 */


#include <cstddef>
#include <utility>

#include "common/AS3Options.hpp"


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
