#include "driver_structure.hpp"



//-----------------------------------------------------------------------------------
// CDriver member functions.
//-----------------------------------------------------------------------------------


CDriver::CDriver
(
 const char *filename
)
 /*
	* Constructor for the driver, which initiates the solver.
	*/
{
	// Initialize a config container.
	mConfigContainer   = std::make_unique<CConfig>(filename);

	// Initialize the geometry container.
	mGeometryContainer = std::make_unique<CGeometry>(mConfigContainer.get());

	// Initialize the temporal container.
	mTemporalContainer = CGenericFactory::CreateTemporalContainer(mConfigContainer.get());

	// Initialize the solver containers.
	mSolverContainer   = CGenericFactory::CreateMultizoneSolverContainer(mConfigContainer.get(), 
			                                                                 mGeometryContainer.get());

	// Initialize the output container.
	mOutputContainer   = std::make_unique<COutput>(mConfigContainer.get(), 
			                                           mGeometryContainer.get());
}

//-----------------------------------------------------------------------------------

CDriver::~CDriver
(
 void
)
 /*
	* Destructor, which cleans up after the driver class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CDriver::InitializeData
(
 void
)
 /*
	* Function that initializes the data for the simulation. 
	*/
{
	// Instantiate an initial condition container. Note, if we are interested in book-keeping 
	// its parameters (e.g. for post-processing or time-dependent initializations), then we 
	// can move this object to a member variable in CDriver.
	std::unique_ptr<IInitialCondition> initial_container;


	// Import an AS3-type grid.
	NImportFile::ImportAS3Grid(mConfigContainer.get(), 
			                       mGeometryContainer.get());

	// Initialize the physical elements, per solver.
	for(size_t iSolver=0; iSolver<mSolverContainer.size(); iSolver++)
	{
		mSolverContainer[iSolver]->InitPhysicalElements(mConfigContainer.get(),
				                                            mGeometryContainer.get());
	}


	// Instantiate the specified initial condition. 
	initial_container = CGenericFactory::CreateInitialConditionContainer(mConfigContainer.get());
	
	// Initialize the solution in each solver/zone.
	for(unsigned short iSolver=0; iSolver<mSolverContainer.size(); iSolver++)
	{
		initial_container->InitializeSolution(mConfigContainer.get(), 
				                                  mGeometryContainer->GetZoneGeometry(iSolver),
																					mSolverContainer[iSolver].get());
	}


	// DEBUGGING
	mOutputContainer->WriteVisualFile(mConfigContainer.get(), 
			                              mGeometryContainer.get(),
																		mSolverContainer);

}







