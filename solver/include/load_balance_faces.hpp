#pragma once

#include "option_structure.hpp"
#include "flattened_index_structure.hpp"


class CLoadBalancedFacePermutation
{
  private:
    
    struct CTypeFaceWeights
    {
      float mWeightInternalFace  = 1.0f;
      float mWeightBoundaryFace  = 1.0f;
      float mWeightInterfaceFace = 1.0f;
    };

  public:
    CLoadBalancedFacePermutation(void) = default;

    CLoadBalancedFacePermutation(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                 size_t                                  nChunk,
                                 EFaceLoadBalanceStrategy                strategy);

    const as3vector1d<size_t> &GetPermutation(void) const { return mPermutation; }
    size_t GetPermutationIndex(size_t i) const
    {
#if DEBUG
      if( i >= mPermutation.size() ) ERROR("Index is out of range.");
#endif
      return mPermutation[i];
    }
    
    size_t GetnFace(void) const { return mPermutation.size(); }

   private:
    CTypeFaceWeights    mWeightsFace;
    as3vector1d<size_t> mPermutation;

    void InitializeFaceWeights(void);
    void InitializeLoadBalancing(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                 size_t                                  nChunk,
                                 EFaceLoadBalanceStrategy                strategy);

    void BuildGreedyPermutation(const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
                                size_t                                  nChunk);

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
