#pragma once

/*!
 * @file
 * @brief Cost-weighted face permutations for static OpenMP work distribution.
 */


#include <cstddef>

#include "option_structure.hpp"
#include "flattened_index_structure.hpp"


/*!
 * @brief A cost-weighted permutation of flattened faces for static worksharing.
 *
 * Construct outside a parallel region. The permutation is reusable while the
 * face collection, actual team size, and schedule(static) ownership stay fixed.
 * Rebuild it if those change; lookup requires no additional synchronization.
 */
class CLoadBalancedFacePermutation
{
  private:
    

    /*!
     * @brief Relative cost estimates used to distribute the three face categories.
     */
    struct CTypeFaceWeights
    {
      float mWeightInternalFace  = 2.0f;  ///< Relative cost assigned to an interior face.
      float mWeightBoundaryFace  = 1.0f;  ///< Relative cost assigned to a physical boundary face.
      float mWeightInterfaceFace = 3.0f;  ///< Relative cost assigned to an interface face.
    };

  public:

    /*!
     * @brief Default constructor that initializes an empty CLoadBalancedFacePermutation.
     */
    CLoadBalancedFacePermutation(void) = default;

    /*!
     * @brief Constructor that initializes CLoadBalancedFacePermutation from the specified data.
     *
     * @param[in] flattened_faces face addresses before applying the load-balancing permutation.
     * @param[in] nChunk requested OpenMP team size used to probe static ownership.
     * @param[in] strategy face load-balancing strategy.
     */
    CLoadBalancedFacePermutation(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                 std::size_t                             nChunk,
                                 EFaceLoadBalanceStrategy                strategy);

    /*!
     * @brief Getter function which returns the original face indices in permuted traversal order.
     *
     * @return original face indices in permuted traversal order.
     */
    const as3vector1d<std::size_t> &GetPermutation(void) const { return mPermutation; }

    /*!
     * @brief Getter function which returns the original flattened face index at the requested permuted
     *   position.
     *
     * @param[in] i zero-based traversal index.
     *
     * @return original flattened face index at the requested permuted position.
     */
    std::size_t GetPermutationIndex(std::size_t i) const
    {
#if DEBUG
      if( i >= mPermutation.size() ) ERROR("Index is out of range.");
#endif
      return mPermutation[i];
    }
    

    /*!
     * @brief Getter function which returns the number of faces.
     *
     * @return number of faces.
     */
    std::size_t GetnFace(void) const { return mPermutation.size(); }

    /*!
     * @brief Actual OpenMP team size observed when constructing the permutation.
     */
    std::size_t GetPlannedTeamSize(void) const { return mPlannedTeamSize; }

    /*!
     * @brief Strategy used to construct this permutation.
     */
    EFaceLoadBalanceStrategy GetStrategy(void) const { return mStrategy; }

    /*!
     * @brief Display name of the strategy used to construct this permutation.
     */
    const char* GetStrategyName(void) const;

    /*!
     * @brief Relative cost of a face category, shared by balancing and reporting.
     */
    float GetFaceWeight(ETypeFaceGeometry face_type) const
    {
      switch( face_type )
      {
        case( ETypeFaceGeometry::INTERNAL ):  { return mWeightsFace.mWeightInternalFace;  }
        case( ETypeFaceGeometry::BOUNDARY ):  { return mWeightsFace.mWeightBoundaryFace;  }
        case( ETypeFaceGeometry::INTERFACE ): { return mWeightsFace.mWeightInterfaceFace; }
        default: ERROR("Unknown face weight.");
      }
    }

   private:
    CTypeFaceWeights         mWeightsFace;  ///< Relative cost estimates for each face category.
    as3vector1d<std::size_t> mPermutation;  ///< Original flattened face index for each position in the permuted traversal.
    std::size_t mPlannedTeamSize = 1;       ///< Actual team size used to probe static slot ownership.
    EFaceLoadBalanceStrategy mStrategy = EFaceLoadBalanceStrategy::GREEDY; ///< Strategy used to build this permutation.

    /*!
     * @brief Function that sets relative cost estimates for each face category.
     */
    void InitializeFaceWeights(void);

    /*!
     * @brief Function that validates the requested distribution and constructs its face permutation.
     *
     * @param[in] flattened_faces face addresses before applying the load-balancing permutation.
     * @param[in] nChunk number of chunks used for static work distribution.
     * @param[in] strategy face load-balancing strategy.
     */
    void InitializeLoadBalancing(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                 std::size_t                             nChunk,
                                 EFaceLoadBalanceStrategy                strategy);

    /*!
     * @brief Function that swaps faces between actual static slots until no single swap reduces imbalance.
     *
     * @param[in] flattened_faces face addresses before applying the load-balancing permutation.
     * @param[in] nChunk requested OpenMP team size used to probe static ownership.
     */
    void BuildGreedyPermutation(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                std::size_t                             nChunk);

};
