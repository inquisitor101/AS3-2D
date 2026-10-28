#include "factory_structure.hpp"
#include "compile_time_dimensions.hpp"

#include <cstddef>
#include <memory>
#include <string>


//-----------------------------------------------------------------------------------
// CGenericFactory member functions.
//-----------------------------------------------------------------------------------


std::unique_ptr<CMonitorData> 
CGenericFactory::CreateMonitoringContainer
(
 const CConfig *config_container
)
 /*
	* Function that constructs the simulation monitoring data container.
	*/
{
	// For now, simply return the only monitoring data container.
	return std::make_unique<CMonitorData>(config_container); 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<ITemporal> 
CGenericFactory::CreateTemporalContainer
(
 const CConfig *config_container
)
 /*
	* Function that constructs the configured temporal method.
	*/
{
	// Check what type of container is specified.
	switch( config_container->GetTemporalScheme() )
	{
		case(ETemporalScheme::SSPRK3):
		{
			return std::make_unique<CSSPRK3Temporal>(config_container);
			break;
		}

		default: ERROR("Unknown type of temporal container.");
	}	

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<IFileVTK> 
CGenericFactory::CreateVTKContainer
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that constructs the configured visualization writer.
	*/
{
	// Check the type of output visualization format.
	switch( config_container->GetOutputVisFormat() )
	{
		// Legacy VTK in binary.
		case( EVisualFormat::VTK_LEGACY_BINARY ):
		{
			return std::make_unique<CLegacyBinaryVTK>(config_container, multizone_geometry_container);
			break;
		}

		default: ERROR("Unknown output visualization format.");
	}

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<IInitialCondition> 
CGenericFactory::CreateInitialConditionContainer
(
 const CConfig *config_container
)
 /*
	* Function that constructs the configured initial-condition evaluator.
	*/
{
	// Check what type of container is specified.
	switch( config_container->GetTypeIC() )
	{
		case(ETypeIC::GAUSSIAN_PRESSURE):
		{
			return std::make_unique<CGaussianPressureIC>(config_container);
			break;
		}

		case(ETypeIC::ISENTROPIC_VORTEX):
		{
			return std::make_unique<CIsentropicVortexIC>(config_container);
			break;
		}

		default: ERROR("Unknown type of initial condition container.");
	}	

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<IRiemannSolver> 
CGenericFactory::CreateRiemannSolverContainer
(
 const CConfig      *config_container,
 ETypeRiemannSolver  riemann
)
 /*
	* Function that constructs the selected numerical surface-flux solver.
	*/
{
	// Check what type of container is specified.
	switch( riemann )
	{
		case(ETypeRiemannSolver::ROE):
		{
			return std::make_unique<CRoeRiemannSolver>(config_container);
			break;
		}

		default: ERROR("Unknown type of Riemann solver container.");
	}	

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<CStandardElement> 
CGenericFactory::CreateStandardElement
(
 const CConfig             *config_container,
 const CSinglezoneGeometry *singlezone_geometry_container
)
 /*
	* Function that constructs the reference element for the zone geometry order and integration rule.
	*/
{
	// Deduce the polynomial order of the solution and type of DOFs. 
	const auto nPolySolution            = singlezone_geometry_container->GetnPolyGrid();
	const auto typeSolutionDistribution = singlezone_geometry_container->GetTypeNodalDistribution(); 

  // Extract the bonus quadrature factor in 1D.
  const auto nBonusQuadratureFactor = config_container->GetBonusQuadratureFactor();

  // Extract the quadrature method.
  const auto typeQuadrature = config_container->GetTypeQuadrature();

  // Extract the type of DG method specified.
  const auto typeDGMethod = config_container->GetTypeDGMethod();

	// Return the correct instantiation of this object.
	return std::make_unique<CStandardElement>(typeSolutionDistribution, 
                                            nPolySolution, 
                                            nBonusQuadratureFactor,
                                            typeQuadrature,
                                            typeDGMethod);
}

//-----------------------------------------------------------------------------------

std::unique_ptr<CPhysicalElement> 
CGenericFactory::CreatePhysicalElement
(
 const CConfig    *config_container,
 CStandardElement *standard_element,
 ITensorProduct   *tensor_container,
 CElementGeometry *element_geometry,
 std::size_t       iZone,
 std::size_t       nVar
)
 /*
	* Function that constructs the physical state and metric storage for one element.
	*/
{
	return std::make_unique<CPhysicalElement>(config_container, 
			                                      standard_element,
																						tensor_container,
																						element_geometry, 
																						iZone, nVar);
}

//-----------------------------------------------------------------------------------

std::unique_ptr<ITensorProduct> 
CGenericFactory::CreateTensorContainer
(
 CStandardElement *standard_element,
 std::size_t    nVar
)
 /*
	* Function that constructs a supported tensor specialization for the reference-element dimensions.
	*/
{
	// Extract the number of solution points (k) and integration points (m) in 1D.
	const std::size_t k = standard_element->GetnSol1D();
	const std::size_t m = standard_element->GetnInt1D();
  const std::size_t n = nVar;

	// Macro which helps in readability for the compile-time specialized tensor classes.
#define SPECIALIZED_TENSOR(K,M,N) if( (K==k) && (M==m) && (N==n) ) \
	return std::make_unique< CTensorProduct<K,M,N> >(standard_element);

	AS3_SUPPORTED_ELEMENT_SPECIALIZATIONS(SPECIALIZED_TENSOR)

#undef SPECIALIZED_TENSOR

	// If the program made it this far, it means the specified values are not implemented.
	ERROR("Combination of (K,M,N) = " + std::to_string(k) + ", " 
                                    + std::to_string(m) + ", "
                                    + std::to_string(n) + " is not found.");			

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<IBoundary> 
CGenericFactory::CreateBoundaryContainer
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *zone_geometry,
 const ISolver              *zone_solver,
 const CBoundaryFaceBatch   &boundary_batch
)
 /*
	* Function that constructs the configured physical boundary evaluator for one face batch.
	*/
{
	// Determine the boundary condition parameters associated with this marker. 
	const auto* param_bc = config_container->GetBoundaryParamMarker( boundary_batch.GetName() ); 

	// Check if the parameter object is found, else issue an error.
	if( !param_bc ) 
	{
		ERROR("Could not find the boundary parameter for the marker: " 
				  + boundary_batch.GetName());
	}


	// Check what type of container is specified.
	switch( param_bc->GetTypeBC() )
	{
    case( ETypeBC::SYMMETRY ):
    {
      return std::make_unique<CSymmetryBoundary>(config_container,
                                                 zone_geometry,
                                                 zone_solver,
                                                 boundary_batch,
                                                 param_bc);
      break;
    }

    case( ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN ):
    {
      return std::make_unique<CRiemannStaticMassInletBoundary>(config_container,
                                                               zone_geometry,
                                                               zone_solver,
                                                               boundary_batch,
                                                               param_bc);
      break;
    }
    
    case( ETypeBC::PRESSURE_OUTLET_RIEMANN ):
    {
      return std::make_unique<CRiemannPressureOutletBoundary>(config_container,
                                                              zone_geometry,
                                                              zone_solver,
                                                              boundary_batch,
                                                              param_bc);
      break;
    }

    default: ERROR("Unknown type of boundary container.");
	}

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<ISolver> 
CGenericFactory::CreateSolverContainer
(
 const CConfig             *config_container,
 const CSinglezoneGeometry *singlezone_geometry
)
 /*
	* Function that constructs the configured zone solver.
	*/
{
  // Extract the zone ID.
  const auto iZone = singlezone_geometry->GetZoneID();

	// Check what type of container is specified.
	switch( config_container->GetTypeSolver(iZone) )
	{
		case(ETypeSolver::EE):
		{
			return std::make_unique<CEESolver>(config_container, singlezone_geometry);
			break;
		}

		default: ERROR("Unknown type of solver container.");
	}	

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<IInterface>
CGenericFactory::CreateInterfaceContainer
(
 const CConfig               *config_container,
 CMultizoneSolver            *multizone_solver_container,
 const CInterfaceFaceBatch   &interface_batch
)
 /*
	* Function that constructs the interface evaluator for a paired-zone face batch.
	*/
{
	// Extract the zone ID of the faces sharing this interface.
	const std::size_t iZone = interface_batch.GetOwnerZoneIndex();
	const std::size_t jZone = interface_batch.GetNeighborZoneIndex();

	// Check what type of solver we have in the iZone.
	switch( multizone_solver_container->GetSinglezoneSolver(iZone)->GetTypeSolver() )
	{
		case(ETypeSolver::EE):
		{
			// Check the type of solver in the jZone too.
			switch( multizone_solver_container->GetSinglezoneSolver(jZone)->GetTypeSolver() )
			{
				// This is a EE-EE interface.
				case(ETypeSolver::EE):
				{
					return std::make_unique<CEEInterface>(config_container, 
																								multizone_solver_container,
                                                interface_batch);
					break;
				}

				default: ERROR("Unknown type of solver container in zone: " + std::to_string(jZone));
			}

			break;
		}

		default: ERROR("Unknown type of solver container in zone: " + std::to_string(iZone));
	}

	// To avoid a compiler warning.
	return nullptr; 
}

