#include "solver/CMultizoneSolver.hpp"

#include "boundary/IBoundary.hpp"
#include "CGenericFactory.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "interface/IInterface.hpp"
#include "solver/ISinglezoneSolver.hpp"

#include <cstddef>


//-----------------------------------------------------------------------------------
// CMultizoneSolver member functions.
//-----------------------------------------------------------------------------------

CMultizoneSolver::CMultizoneSolver
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Constructor that creates the solver for each grid zone.
  */
{
  // Deduce the number of zones.
  mNZone = multizone_geometry_container->GetnZone();

	// Allocate the necessary number of solvers.
	mMultizoneSolverContainer.reserve( mNZone );

	// Initialize each solver.
	for(std::size_t iZone=0; iZone<mNZone; iZone++)
	{
    // Consistency check.
    if( iZone != multizone_geometry_container->GetSinglezoneGeometry(iZone)->GetZoneID() )
    {
      ERROR("Zone indices do not match.");
    }

    // Instantiate the current solver in this zone.
		mMultizoneSolverContainer.emplace_back
    (
     CGenericFactory::CreateSolverContainer(config_container, multizone_geometry_container->GetSinglezoneGeometry(iZone) )
    );
	}
}

//-----------------------------------------------------------------------------------

CMultizoneSolver::~CMultizoneSolver(void) = default;

//-----------------------------------------------------------------------------------

void CMultizoneSolver::InitializeInterfaces
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Function that creates interface evaluators and their direction-local batch mappings.
  */
{
  // Extract the interfaces.
  const auto& idir_interfaces = multizone_geometry_container->GetInterfaceBatchesIFace();
  const auto& jdir_interfaces = multizone_geometry_container->GetInterfaceBatchesJFace();

  // Reserve the required memory for the indicial mapping of the interfaces.
  mIndexInterfaceIFace.reserve( idir_interfaces.size() );
  mIndexInterfaceJFace.reserve( jdir_interfaces.size() );
  
  // Reserve the required number of (unique) interface containers.
  mMultizoneInterfaceContainer.reserve( idir_interfaces.size() + jdir_interfaces.size() );


  // Loop over each interface which contains the ith face as a IFace.
  for(const auto& batch : idir_interfaces )
  {
    // Global-book-keeping index.
    const std::size_t iInterface = mMultizoneInterfaceContainer.size();

    // Initialize each interface object.
    mMultizoneInterfaceContainer.emplace_back
    (
     CGenericFactory::CreateInterfaceContainer(config_container,
                                               this, batch)
    );

    // Book-keep the index map.
    mIndexInterfaceIFace.push_back( iInterface );
  }


  // Loop over each interface which contains the ith face as a IFace.
  for(const auto& batch : jdir_interfaces )
  {
    // Global-book-keeping index.
    const std::size_t iInterface = mMultizoneInterfaceContainer.size();

    // Initialize each interface object.
    mMultizoneInterfaceContainer.emplace_back
    (
     CGenericFactory::CreateInterfaceContainer(config_container,
                                               this, batch)
    );
  
    // Book-keep the index map.
    mIndexInterfaceJFace.push_back( iInterface );
  }
}

//-----------------------------------------------------------------------------------

void CMultizoneSolver::InitializeBoundaries
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
  * Function that creates boundary evaluators and their direction-local batch mappings.
  */
{
  // Extract the boundaries.
  const auto& idir_boundaries = multizone_geometry_container->GetBoundaryBatchesIFace();
  const auto& jdir_boundaries = multizone_geometry_container->GetBoundaryBatchesJFace();

  // Reserve the required memory for the indicial mapping of the boundaries.
  mIndexBoundaryIFace.reserve( idir_boundaries.size() );
  mIndexBoundaryJFace.reserve( jdir_boundaries.size() );
  
  // Reserve the required number of boundary containers.
  mMultizoneBoundaryContainer.reserve( idir_boundaries.size() + jdir_boundaries.size() );

  // Loop over each boundary which contains with an IFace.
  for(const auto& batch : idir_boundaries )
  {
    // Global-book-keeping index.
    const std::size_t iBoundary = mMultizoneBoundaryContainer.size();

    // Extract the current zone index.
    const std::size_t iZone = batch.GetZoneIndex();

    // Extract the current zone solver and geometry.
    const auto* zone_solver   = GetSinglezoneSolver(iZone);
    const auto* zone_geometry = multizone_geometry_container->GetSinglezoneGeometry(iZone);

    // Initialize each boundary object.
    mMultizoneBoundaryContainer.emplace_back
    (
     CGenericFactory::CreateBoundaryContainer(config_container,
                                              zone_geometry,
                                              zone_solver, 
                                              batch)
    );

    // Book-keep the index map.
    mIndexBoundaryIFace.push_back( iBoundary );
  }


  // Loop over each boundary which contains with an JFace.
  for(const auto& batch : jdir_boundaries )
  {
    // Global-book-keeping index.
    const std::size_t iBoundary = mMultizoneBoundaryContainer.size();

    // Extract the current zone index.
    const std::size_t iZone = batch.GetZoneIndex();

    // Extract the current zone solver and geometry.
    const auto* zone_solver   = GetSinglezoneSolver(iZone);
    const auto* zone_geometry = multizone_geometry_container->GetSinglezoneGeometry(iZone);

    // Initialize each boundary object.
    mMultizoneBoundaryContainer.emplace_back
    (
     CGenericFactory::CreateBoundaryContainer(config_container,
                                              zone_geometry,
                                              zone_solver, 
                                              batch)
    );

    // Book-keep the index map.
    mIndexBoundaryJFace.push_back( iBoundary );
  }
}
