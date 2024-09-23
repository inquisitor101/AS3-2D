#include "factory_structure.hpp"


//-----------------------------------------------------------------------------------
// CGenericFactory member functions.
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
 CConfig   *config_container,
 CGeometry *geometry_container
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
 CConfig           *config_container,
 ETypeRiemannSolver riemann
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
 CConfig       *config_container,
 unsigned short iZone
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
 CStandardElement *standard_element
)
 /*
	* Function that creates a specialized instance of a templated tensor container.
	*/
{
	// Extract the number of solution points (k) and integration points (m) in 1D.
	const size_t k = standard_element->GetnSol1D();
	const size_t m = standard_element->GetnInt1D();

	// Macro which helps in readability for the compile-time specialized tensor classes.
#define SPECIALIZED_TENSOR(K,M) if( (K==k) && (M==m) ) \
	return std::make_unique< CTensorProduct<K,M> >(standard_element);

	SPECIALIZED_TENSOR(2,3);

	SPECIALIZED_TENSOR(3,3);
	SPECIALIZED_TENSOR(3,4);

	SPECIALIZED_TENSOR(4,4);
	SPECIALIZED_TENSOR(4,5);

	SPECIALIZED_TENSOR(5,5);
	SPECIALIZED_TENSOR(5,7);

	SPECIALIZED_TENSOR(6,6);
	SPECIALIZED_TENSOR(6,8);

	SPECIALIZED_TENSOR(7, 7);
	SPECIALIZED_TENSOR(7,10);

	SPECIALIZED_TENSOR(8, 8);
	SPECIALIZED_TENSOR(8,11);

	SPECIALIZED_TENSOR(9, 9);
	SPECIALIZED_TENSOR(9,13);
	
	SPECIALIZED_TENSOR(10,10);
	SPECIALIZED_TENSOR(10,14);

	// If the program made it this far, it means the specified values are not implemented.
	ERROR("Combination of (K,M) = " + std::to_string(k) + ", " + std::to_string(m) + " is not found.");			

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
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
 /*
	* Function that creates a specialized instance of a solver container.
	*/
{
	// Check what type of container is specified.
	switch( config_container->GetTypeSolver(iZone) )
	{
		case(ETypeSolver::EE):
		{
			return std::make_unique<CEESolver>(config_container, geometry_container, iZone);
			break;
		}

		default: ERROR("Unknown type of solver container.");
	}	

	// To avoid a compiler warning.
	return nullptr; 
}

//-----------------------------------------------------------------------------------

as3vector1d<std::unique_ptr<ISolver>>
CGenericFactory::CreateMultizoneSolverContainer
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that creates a vector of specialized instances of solver containers.
	*/
{
	// Allocate the necessary number of solvers.
	as3vector1d<std::unique_ptr<ISolver>> solver_container( config_container->GetnZone() );

	// Initialize each solver.
	for(unsigned short iZone=0; iZone<solver_container.size(); iZone++)
	{
		solver_container[iZone] = CreateSolverContainer(config_container, geometry_container, iZone);
	}

	// Return the vector of solver containers.
	return solver_container; 
}

//-----------------------------------------------------------------------------------

std::unique_ptr<IInterface>
CGenericFactory::CreateInterfaceContainer
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 CInterfaceParamMarker                 *param_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container 
)
 /*
	* Function that creates a specialized instance of an interface boundary container.
	*/
{
	// Temporary lambda to search for the zone of a given marker name.
	auto lFindMarker = [=](std::string &name) -> CMarker*
	{
		for( auto& zone: geometry_container->GetZoneGeometry() )
		{
			for( auto& marker: zone->GetMarker() )
			{
				if( marker->GetNameMarker() == name )
				{
					return marker.get();
				}
			}
		}
		
		// The program should not reach here, otherwise we have an error.
		ERROR("Could not find the interface marker.");

		// To avoid compiler problems, return something.
		return nullptr;
	};

	// Get a pointer to the two markers forming this interface.
	const CMarker *imarker_container = lFindMarker(param_container->mName);
	const CMarker *jmarker_container = lFindMarker(param_container->mNameMatching);

	// Extract the zone ID of these markers.
	const unsigned short iZone = imarker_container->GetZoneID();
	const unsigned short jZone = jmarker_container->GetZoneID();

	// Check what type of solver we have in the iZone.
	switch( solver_container[iZone]->GetTypeSolver() )
	{
		case(ETypeSolver::EE):
		{
			// Check the type of solver in the jZone too.
			switch( solver_container[jZone]->GetTypeSolver() )
			{
				// This is a EE-EE interface.
				case(ETypeSolver::EE):
				{
					return std::make_unique<CEEInterface>(config_container, 
							                                  geometry_container, 
																								param_container,
																								imarker_container, 
																								jmarker_container,
																								solver_container);
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



