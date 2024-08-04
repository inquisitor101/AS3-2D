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
	
	// Instantiate physical element objects.
	for(size_t i=0; i<mPhysicalElementContainer.size(); i++)
	{
		mPhysicalElementContainer[i] = CGenericFactory::CreatePhysicalElement(config_container,
				                                                                  mStandardElementContainer.get(),
																																					mTensorProductContainer.get(),
				                                                                  zone->GetElementGeometry(i),
																																					mZoneID, mNVar);
	}
}




