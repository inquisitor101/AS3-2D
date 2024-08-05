#include "solver_structure.hpp"


//-----------------------------------------------------------------------------------
// ISolver member functions.
//-----------------------------------------------------------------------------------


ISolver::ISolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
	:
		mZoneID(iZone)
 /*
	* Constructor for the interface solver class.
	*/
{
	// Instantiate a standard element object.
	mStandardElementContainer = CGenericFactory::CreateStandardElement(config_container, iZone);

	// Instantiate a tensor-product object.
	mTensorProductContainer = CGenericFactory::CreateTensorContainer(mStandardElementContainer.get());

	// Instantiate a spatial object.
	mSpatialContainer = CGenericFactory::CreateSpatialContainer(config_container, geometry_container, iZone);

	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(iZone);

	// Ensure the zone matches this one.
	if( iZone != zone->GetZoneID() ) ERROR("Geometry zone does not match Solver zone.");
}

//-----------------------------------------------------------------------------------

ISolver::~ISolver
(
 void
)
 /*
	* Destructor, which cleans up after the interface solver class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CEESolver member functions.
//-----------------------------------------------------------------------------------


CEESolver::CEESolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
	:
		ISolver(config_container, geometry_container, iZone)
 /*
	* Constructor for the (non-linear) Euler equations class.
	*/
{

}

//-----------------------------------------------------------------------------------

CEESolver::~CEESolver
(
 void
)
 /*
	* Destructor, which cleans up after the Euler equations class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEESolver::InitPhysicalElements
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the physical elements. 
	*/
{
	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(mZoneID);

	// Allocate memory for the physical elements in this solver.
	mPhysicalElementContainer.resize( zone->GetnElem() );
	
	// Instantiate physical element objects. The reason this needs a factory is because we might have 
	// a specialized solver (e.g. EE), but with a buffer zone, such as a PML. Hence, the elements
	// require additional variables. For now, avoid an interface class, as this object is the most 
	// performance-reliant data class in the code. Perhaps it is better to keep a single class, but 
	// instantiate it differently, depending on the type of solver and buffer layer specification.
	for(size_t i=0; i<mPhysicalElementContainer.size(); i++)
	{
		mPhysicalElementContainer[i] = CGenericFactory::CreatePhysicalElement(config_container,
				                                                                  mStandardElementContainer.get(),
																																					mTensorProductContainer.get(),
				                                                                  zone->GetElementGeometry(i),
																																					mZoneID, mNVar);
	}
}





