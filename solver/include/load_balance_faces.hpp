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
     * @param[in] nChunk number of chunks used for static work distribution.
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

   private:
    CTypeFaceWeights         mWeightsFace;  ///< Relative cost estimates for each face category.
    as3vector1d<std::size_t> mPermutation;  ///< Original flattened face index for each position in the permuted traversal.

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
     * @brief Function that assigns costly faces first to the least-loaded chunk and forms a static traversal.
     *
     * @param[in] flattened_faces face addresses before applying the load-balancing permutation.
     * @param[in] nChunk number of chunks used for static work distribution.
     */
    void BuildGreedyPermutation(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                std::size_t                             nChunk);

    /*!
     * @brief Getter function which returns the relative cost associated with a face category.
     *
     * @param[in] face_type interior, boundary or interface face category.
     *
     * @return relative cost associated with a face category.
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

};
