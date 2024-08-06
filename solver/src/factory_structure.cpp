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

std::unique_ptr<ISpatial> 
CGenericFactory::CreateSpatialContainer
(
 CConfig        *config_container,
 CGeometry      *geometry_container,
 unsigned short  iZone
)
 /*
	* Function that creates a specialized instance of a spatial container.
	*/
{
	// Check what type of container is specified.
	switch( config_container->GetTypeSolver(iZone) )
	{
		case(ETypeSolver::EE):
		{
			switch( config_container->GetTypeBufferLayer(iZone) )
			{
				case(ETypeBufferLayer::NONE): 
				{
					return std::make_unique<CEESpatial>(config_container, geometry_container);
					break;
				}
				default: ERROR("Unknown type of buffer layer.");
			}

			break;
		} // End: EE

		default: ERROR("Unknown type of solver.");
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
	return std::make_unique<CStandardElement>(config_container, iZone);
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
	* The instantiation is based on the input standard_element, namely <K,M>: 
	*   K is equivalent to nSol1D.
	*   M is equivalent to nInt1D.
	*/
{
	// While we can create a static variable, it doesn't make sense to have all the 
	// different objects instantiated over the entire lifetime of the simulation. In 
	// principle, we can use shared_ptr with static members, but creating these objects
	// and then deleting all except the one corresponding to <m,k> is fine as a 
	// preprocessing step.
	std::vector<std::unique_ptr<ITensorProduct>> mTensor;

	// Lambda that accumulates the specialized templated classes, based on dimensions M and K.
	auto f = [&mTensor]<size_t kk, size_t mm>() -> void
	{
		mTensor.push_back( std::make_unique<CTensorProduct<kk,mm>>() );
	};

	// Macro to help write cleaner code, since lambda-template syntax is (very) ugly..
#define ADD_TENSOR(KK,MM) f.template operator()<KK,MM>()


	// K = 2.
	ADD_TENSOR(2,2); // M = 2.
	
	// K = 3.
	ADD_TENSOR(3,3); // M = 3.
	ADD_TENSOR(3,4); // M = 4.
	
	// K = 4.
	ADD_TENSOR(4,4); // M = 4.
	ADD_TENSOR(4,5); // M = 5.
	
	// K = 5.
	ADD_TENSOR(5,5); // M = 5.
	ADD_TENSOR(5,7); // M = 7.
	
	// K = 6.
	ADD_TENSOR(6,6); // M = 6.
	ADD_TENSOR(6,8); // M = 8.
	
	// K = 7.
	ADD_TENSOR(7, 7); // M = 7.
	ADD_TENSOR(7,10); // M = 10.
	
	// K = 8.
	ADD_TENSOR(8, 8); // M = 8.
	ADD_TENSOR(8,11); // M = 11.
	
	// K = 9.
	ADD_TENSOR(9, 9); // M = 9.
	ADD_TENSOR(9,13); // M = 13.
	
	// K = 10.
	ADD_TENSOR(10,10); // M = 10.
	ADD_TENSOR(10,14); // M = 14.


	// Check if value exists and move its ownership, otherwise break the code.
	const size_t k = standard_element->GetnSol1D();
	const size_t m = standard_element->GetnInt1D();
	for(auto& p: mTensor)		
	{
		if( (p->mK==k) && (p->mM==m) )
		{
			// Assign its standard element pointer.
			p->SetStandardElement(standard_element);
			return std::move(p);
		}
	}

	// If the code made it this far, it means the values have not been implemented (yet). 
	ERROR("Combination of (K,M) = " + std::to_string(k) + ", " + std::to_string(m) + " is not found.");

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
	as3vector1d<std::unique_ptr<ISolver>> solver;

	// Check what type of container is specified.
	for(unsigned short iZone=0; iZone<config_container->GetnZone(); iZone++)
	{
		solver.push_back( CreateSolverContainer(config_container, geometry_container, iZone) );
	}

	// To avoid a compiler warning.
	return solver; 
}
