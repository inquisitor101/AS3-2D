#include "load_balance_faces.hpp"

#include <cstddef>
#include <iostream>
#include <ostream>
#include <array>
#include <limits>
#include <numeric>
#include <set>
#include <algorithm>


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

const char* CLoadBalancedFacePermutation::GetStrategyName(void) const
{
  switch( mStrategy )
  {
    case EFaceLoadBalanceStrategy::GREEDY: return "Greedy";
    default: ERROR("Unknown face load-balancing strategy.");
  }
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
  mStrategy = strategy;
  switch( mStrategy )
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
  * Function that balances actual static slots by exchanging faces between workers.
  */
{
  const std::size_t nFace = flattened_faces.size();

  // Keep the original order as the starting point. Every accepted swap reduces
  // squared load imbalance, and neither affected worker exceeds its old maximum.
  mPermutation.resize(nFace);
  std::iota(mPermutation.begin(), mPermutation.end(), std::size_t{0});

  // Face categories index the slot lists; their costs come from the same model
  // used elsewhere in this class, rather than assuming costs are enum values.
  const std::array<float, 3> weights = {
    GetFaceWeight(ETypeFaceGeometry::BOUNDARY),
    GetFaceWeight(ETypeFaceGeometry::INTERNAL),
    GetFaceWeight(ETypeFaceGeometry::INTERFACE)
  };
  const auto category = [](ETypeFaceGeometry type) -> std::size_t
  {
    switch( type )
    {
      case ETypeFaceGeometry::BOUNDARY: return 0;
      case ETypeFaceGeometry::INTERNAL: return 1;
      case ETypeFaceGeometry::INTERFACE: return 2;
      default: ERROR("Unknown face type.");
    }
  };

#ifndef HAVE_OPENMP
  nChunk = 1;
#else
  if( nChunk > static_cast<std::size_t>(std::numeric_limits<int>::max()) )
    ERROR("OpenMP thread count exceeds the supported range.");
#endif

  // Probe the actual schedule once. Do not cap the team at nFace: doing so can
  // change which worker owns each slot when there are more workers than faces.
  // Each thread writes only its own slot lists and load.
  using CCategorySlots = std::array<as3vector1d<std::size_t>, 3>;
  as3vector1d<CCategorySlots> slots(nChunk);
  as3vector1d<double> loads(nChunk, 0.0);
  mPlannedTeamSize = 1;
#ifdef HAVE_OPENMP
#pragma omp parallel num_threads(static_cast<int>(nChunk))
#endif
  {
    std::size_t thread = 0;
#ifdef HAVE_OPENMP
    thread = static_cast<std::size_t>(omp_get_thread_num());
#pragma omp single
    mPlannedTeamSize = static_cast<std::size_t>(omp_get_num_threads());
#pragma omp for schedule(static)
#endif
    for(std::size_t slot=0; slot<nFace; ++slot)
    {
      const auto type = category(flattened_faces[slot].mFaceType);
      slots[thread][type].push_back(slot);
      loads[thread] += weights[type];
    }
  }
  slots.resize(mPlannedTeamSize);
  loads.resize(mPlannedTeamSize);

  // For each category, order the workers holding such faces by total load.
  // The best pair for a high/low cost combination is found at the two ends.
  // This avoids searching all faces or all worker pairs for every swap.
  using CLoadEntry = std::pair<double, std::size_t>;
  std::array<std::set<CLoadEntry>, 3> workers;
  const auto updateWorker = [&](std::size_t thread, bool insert)
  {
    for(std::size_t type=0; type<weights.size(); ++type)
      if( !slots[thread][type].empty() )
      {
        const CLoadEntry entry{loads[thread], thread};
        if( insert ) workers[type].insert(entry);
        else workers[type].erase(entry);
      }
  };
  for(std::size_t thread=0; thread<mPlannedTeamSize; ++thread)
    updateWorker(thread, true);

  for(;;)
  {
    double bestImprovement = 0.0;
    std::size_t heavy = 0, light = 0, expensive = 0, cheap = 0;
    for(std::size_t high=0; high<weights.size(); ++high)
      for(std::size_t low=0; low<weights.size(); ++low)
      {
        const double delta = static_cast<double>(weights[high]) - weights[low];
        if( delta <= 0.0 || workers[high].empty() || workers[low].empty() ) continue;
        const auto& mostLoaded = *workers[high].rbegin();
        const auto& leastLoaded = *workers[low].begin();
        const double gap = mostLoaded.first - leastLoaded.first;
        // Half the decrease in the sum of squared thread loads.
        const double improvement = delta*(gap - delta);
        if( improvement > bestImprovement )
        {
          bestImprovement = improvement;
          heavy = mostLoaded.second;
          light = leastLoaded.second;
          expensive = high;
          cheap = low;
        }
      }
    if( bestImprovement <= 0.0 ) break;

    updateWorker(heavy, false);
    updateWorker(light, false);
    const auto heavySlot = slots[heavy][expensive].back();
    const auto lightSlot = slots[light][cheap].back();
    slots[heavy][expensive].pop_back();
    slots[light][cheap].pop_back();
    std::swap(mPermutation[heavySlot], mPermutation[lightSlot]);
    slots[heavy][cheap].push_back(heavySlot);
    slots[light][expensive].push_back(lightSlot);
    const double delta = static_cast<double>(weights[expensive]) - weights[cheap];
    loads[heavy] -= delta;
    loads[light] += delta;
    updateWorker(heavy, true);
    updateWorker(light, true);
  }
}

//-----------------------------------------------------------------------------------







