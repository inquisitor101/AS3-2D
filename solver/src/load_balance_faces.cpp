#include "load_balance_faces.hpp"

#include <cstddef>
#include <iostream>
#include <ostream>
#include <queue>
#include <algorithm>
#include <functional>


//-----------------------------------------------------------------------------------
// CLoadBalancedFacePermutation member functions.
//-----------------------------------------------------------------------------------

CLoadBalancedFacePermutation::CLoadBalancedFacePermutation
(
 const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
 std::size_t                             nChunk,
 EFaceLoadBalanceStrategy                strategy
)
 /*
  * Constructor that builds a cost-weighted face traversal.
  */
{ 
  InitializeLoadBalancing(flattened_faces, nChunk, strategy);
}

//-----------------------------------------------------------------------------------

void CLoadBalancedFacePermutation::InitializeFaceWeights
(
 void
)
 /*
  * Function that sets relative costs for interior, boundary and interface faces.
  */
{
  // NOTE, for now this is hard-coded, but this can be adjusted based on
  // benchmark runs, such that we obtain a normalized estimate with
  // respect to the internal face cost.
  mWeightsFace.mWeightInternalFace  = 2.0f;
  mWeightsFace.mWeightBoundaryFace  = 1.0f; // this needs to be adjusted based on the BC for optimal results.
  mWeightsFace.mWeightInterfaceFace = 3.0f;
}

//-----------------------------------------------------------------------------------

void CLoadBalancedFacePermutation::InitializeLoadBalancing
(
 const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
 std::size_t                             nChunk,
 EFaceLoadBalanceStrategy                strategy
)
 /*
  * Function that validates the requested distribution and builds the face permutation.
  */
{
  // Consistency check.
  if( flattened_faces.size() == 0 ) ERROR("Detected no faces.");
  if( nChunk == 0 ) ERROR("Number of chunks must be greater than zero.");

  // Initialize the face weights, based on the type of faces.
  InitializeFaceWeights();
  
  // Reset and initialize the permutation size.
  mPermutation.clear();

  // Decide which strategy to use.
  switch( strategy )
  {
    case( EFaceLoadBalanceStrategy::GREEDY ):
    {
      BuildGreedyPermutation(flattened_faces, nChunk);
      break;
    }

    default: ERROR("Unknown strategy specified.");
  }
}

//-----------------------------------------------------------------------------------

void CLoadBalancedFacePermutation::BuildGreedyPermutation
(
 const as3vector1d<CFlattenedFaceIndex> &flattened_faces,
 std::size_t                             nChunk
)
 /*
  * Function that assigns costly faces first to the least-loaded traversal chunk.
  */
{
  // Storage for faces and their weights.
  struct CWeightedFace
  {
    std::size_t mIndex;
    float       mWeight;
  };

  // Deduce the number of faces.
  const std::size_t nFace = flattened_faces.size();

  // Issue warning, in case of inefficient implementation.
  if( nChunk > nFace )
  {
    WARNING("Inefficient implementation: nChunk > nFace.");
    nChunk = nFace; // Maybe issue an error instead?
  }

  // Create and reserve memory for the weighted faces.
  as3vector1d<CWeightedFace> weighted_faces;
  weighted_faces.reserve( nFace );

  // Loop over each face and assign its weight.
  for(std::size_t iFace=0; iFace<nFace; iFace++)
  {
    const auto weight = GetFaceWeight( flattened_faces[iFace].mFaceType );
    weighted_faces.emplace_back( CWeightedFace{iFace, weight} );
  } 

  // Sort the most expensive faces first. This is the Longest Processing Time (LPT) ordering 
  // and improves the quality of the greedy partitioning.
  std::sort( weighted_faces.begin(), weighted_faces.end(),
             [](const CWeightedFace &lhs, const CWeightedFace &rhs) { return lhs.mWeight > rhs.mWeight; } );


  // Create a temporary object to store the load of each chunk.
  struct CChunkLoad
  {
    std::size_t mIndexChunk;
    float       mWeight;

    bool operator>(const CChunkLoad& other) const
     /*
      * Function that compares chunk weights for the minimum-load priority queue.
      */
    {
      return mWeight > other.mWeight;
    }
  };

  std::priority_queue<CChunkLoad, as3vector1d<CChunkLoad>, std::greater<CChunkLoad>> chunk_queue;

  // Initialize all chunks with zero load.
  for(std::size_t iChunk=0; iChunk<nChunk; iChunk++)
  {
    chunk_queue.push( CChunkLoad{iChunk, 0.0f} );
  }

  // Store the faces belonging to each chunk.
  as3vector2d<std::size_t> chunk_faces(nChunk);
  as3vector1d<float>  chunk_weights(nChunk, 0.0f);


  // Proceed with the greedy assignment. The next most expensive face is assigned 
  // to the currently least-loaded chunk.
  for(const auto& face : weighted_faces)
  {
    auto chunk = chunk_queue.top();
    chunk_queue.pop();

    chunk_faces[chunk.mIndexChunk].push_back( face.mIndex );
    chunk_weights[chunk.mIndexChunk] += face.mWeight;
    chunk.mWeight += face.mWeight;
    chunk_queue.push(chunk);
  }

  // Finally, construct the final permutation. Each chunk is kept contiguous so that the 
  // chunks (e.g. threads in OpenMP schedule(static)) assigns one balanced region to each thread.
  mPermutation.clear();
  mPermutation.reserve( nFace );

  for(std::size_t iChunk=0; iChunk<nChunk; iChunk++)
  {
    mPermutation.insert( mPermutation.end(), 
                         chunk_faces[iChunk].begin(),
                         chunk_faces[iChunk].end() );
  }

  // Display progress.
  std::cout << "\nReporting load balancing progress:\n";
  for(std::size_t iChunk = 0; iChunk<nChunk; iChunk++)
  {
    std::cout
        << "  load[" << iChunk << "]: "
        << chunk_weights[iChunk]
        << "\n";
  }
  std::cout << std::endl;
}

//-----------------------------------------------------------------------------------







