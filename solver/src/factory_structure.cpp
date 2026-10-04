#include "factory_structure.hpp"


//-----------------------------------------------------------------------------------
// CGenericFactory member functions.
//-----------------------------------------------------------------------------------


std::unique_ptr<CMonitorData> 
CGenericFactory::CreateMonitoringContainer
(
 CConfig *config_container
)
 /*
	* Function that creates a specialized instance of a temporal container.
	*/
{
	// For now, simply return the only monitoring data container.
	return std::make_unique<CMonitorData>(config_container); 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<ITemporal> 
CGenericFactory::CreateTemporalContainer
(
 CConfig *config_container
)
 /*
	* Function that creates a specialized instance of a temporal container.
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
 const CConfig   *config_container,
 const CGeometry *geometry_container
)
 /*
	* Function that creates a specialized instance of a vtk container.
	*/
{
	// Check the type of output visualization format.
	switch( config_container->GetOutputVisFormat() )
	{
		// Legacy VTK in binary.
		case( EVisualFormat::VTK_LEGACY_BINARY ):
		{
			return std::make_unique<CLegacyBinaryVTK>(config_container, geometry_container);
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
 CConfig *config_container
)
 /*
	* Function that creates a specialized instance of an initial condition container.
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
	* Function that creates a specialized instance of a Riemann solver container.
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
 const CConfig  *config_container,
 unsigned short  iZone
)
 /*
	* Function that creates a specialized instance of a standard element container.
	*/
{
	// Deduce the polynomial order of the solution and type of DOFs. 
	auto nPoly   = config_container->GetnPoly(iZone);
	auto typeDOF = config_container->GetTypeDOF(iZone);
	
	// Get the number of (over-)integration points in 1D.
	auto nInt1D  = NPolynomialUtility::IntegrationRule1D(nPoly); 

	// Return the correct instantiation of this object.
	return std::make_unique<CStandardElement>(typeDOF, nPoly, nInt1D);
}

//-----------------------------------------------------------------------------------

std::unique_ptr<CPhysicalElement> 
CGenericFactory::CreatePhysicalElement
(
 CConfig          *config_container,
 CStandardElement *standard_element,
 ITensorProduct   *tensor_container,
 CElementGeometry *element_geometry,
 unsigned short    iZone,
 unsigned short    nVar
)
 /*
	* Function that creates a specialized instance of a physical element container.
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
 unsigned short    nVar
)
 /*
	* Function that creates a specialized instance of a templated tensor container.
	*/
{
	// Extract the number of solution points (k) and integration points (m) in 1D.
	const size_t k = standard_element->GetnSol1D();
	const size_t m = standard_element->GetnInt1D();
  const size_t n = nVar;

	// Macro which helps in readability for the compile-time specialized tensor classes.
#define SPECIALIZED_TENSOR(K,M,N) if( (K==k) && (M==m) && (N==n) ) \
	return std::make_unique< CTensorProduct<K,M,N> >(standard_element);

	SPECIALIZED_TENSOR(2,2,4);
	SPECIALIZED_TENSOR(2,3,4);

	SPECIALIZED_TENSOR(3,3,4);
	SPECIALIZED_TENSOR(3,4,4);

	SPECIALIZED_TENSOR(4,4,4);
	SPECIALIZED_TENSOR(4,5,4);

	SPECIALIZED_TENSOR(5,5,4);
	SPECIALIZED_TENSOR(5,7,4);

	SPECIALIZED_TENSOR(6,6,4);
	SPECIALIZED_TENSOR(6,8,4);

	SPECIALIZED_TENSOR(7, 7,4);
	SPECIALIZED_TENSOR(7,10,4);

	SPECIALIZED_TENSOR(8, 8,4);
	SPECIALIZED_TENSOR(8,11,4);

	SPECIALIZED_TENSOR(9, 9,4);
	SPECIALIZED_TENSOR(9,13,4);
	
	SPECIALIZED_TENSOR(10,10,4);
	SPECIALIZED_TENSOR(10,14,4);

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
 CConfig      *config_container,
 CGeometry    *geometry_container,
 CMarker      *marker_container,
 unsigned int  index
)
 /*
	* Function that creates a specialized instance of a boundary container.
	*/
{
	// Determine the boundary condition associated with this marker. 
	auto* param = config_container->GetBoundaryParamMarker( marker_container->GetNameMarker() ); 

	// Check if the parameter object is found, else issue an error.
	if( !param ) 
	{
		ERROR("Could not find the boundary parameter for the marker: " 
				  + marker_container->GetNameMarker() 
			    + " with element index: " + std::to_string(index) );
	}


	// Check what type of container is specified.
	switch( param->GetTypeBC() )
	{
		
		default: ERROR("Unknown type of boundary container.");
	}

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<ISolver> 
CGenericFactory::CreateSolverContainer
(
 const CConfig       *config_container,
 const CZoneGeometry *zone_geometry
)
 /*
	* Function that creates a specialized instance of a solver container.
	*/
{
  // Extract the zone ID.
  const auto iZone = zone_geometry->GetZoneID();

	// Check what type of container is specified.
	switch( config_container->GetTypeSolver(iZone) )
	{
		case(ETypeSolver::EE):
		{
			return std::make_unique<CEESolver>(config_container, zone_geometry);
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
 const CInterfaceFacesFamily &interface_family
)
 /*
	* Function that creates a specialized instance of an interface boundary container.
	*/
{
	// Extract the zone ID of the faces sharing this interface.
	const unsigned short iZone = interface_family.GetiZone();
	const unsigned short jZone = interface_family.GetjZone();

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
                                                interface_family);
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

