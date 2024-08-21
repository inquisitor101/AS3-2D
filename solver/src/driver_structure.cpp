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
	mConfigContainer    = std::make_unique<CConfig>(filename);

	// Initialize the geometry container.
	mGeometryContainer  = std::make_unique<CGeometry>(mConfigContainer.get());

	// Initialize the specified initial condition. 
	mInitialContainer   = CGenericFactory::CreateInitialConditionContainer(mConfigContainer.get());
	
	// Initialize the temporal container.
	mTemporalContainer  = CGenericFactory::CreateTemporalContainer(mConfigContainer.get());

	// Initialize the solver containers.
	mSolverContainer    = CGenericFactory::CreateMultizoneSolverContainer(mConfigContainer.get(), 
			                                                                  mGeometryContainer.get());

	// Initialize the output container.
	mOutputContainer    = std::make_unique<COutput>(mConfigContainer.get(), 
			                                            mGeometryContainer.get());

	// Initialize the iteration container, must be initialized after the solver container.
	mIterationContainer = std::make_unique<CIteration>(mConfigContainer.get(),
			                                               mSolverContainer); 
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

void CDriver::StartSolver
(
 void
)
 /*
	* Function that runs the entire simulation. 
	*/
{
	// Record start time, based on processor time (not physical time).
	const as3double proc_t0 = as3double( std::clock() )/as3double( CLOCKS_PER_SEC );
	// Record start time, based on the actual physical time.
	const auto phys_t0 = std::chrono::high_resolution_clock::now();


	// Preprocess the data once, prior to starting the simulation. 
	PreProcess();

	// Fire up the actual solver.
	Run();


	// Record the end time used by the processor (not physical time).
	const as3double proc_t1 = as3double( std::clock() )/as3double( CLOCKS_PER_SEC );
	// Record the end time, based on the actual physical time.
	const auto phys_t1 = std::chrono::high_resolution_clock::now();

	// Lapse time used by the entire solver, in processor time (not physical time).
	const as3double lapsedtime_proc = proc_t1 - proc_t0;
	// Lapse time used by the entire solver, in physical time.
	const auto lapsedtime_phys = std::chrono::duration_cast<std::chrono::milliseconds>(phys_t1-phys_t0);

	// Report lapsed time.
	std::cout << "\n% % % % % % % % % % % % % % % % % % % % % % %" << std::endl;
	std::cout << std::scientific << std::setprecision(10) << std::setw(10)
	          << "lapsed (proc) time [sec]: " << lapsedtime_proc                << "  %" << "\n"
	          << "lapsed (phys) time [sec]: " << lapsedtime_phys.count()*1.0e-3 << "  %" << std::endl;
}

//-----------------------------------------------------------------------------------

void CDriver::PreProcess
(
 void
)
 /*
	* Function that preprocesses the data, prior to starting the simulation.
	*/
{
	// First, initialize the solution and physical elements.
	InitializeData();
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
	// Import an AS3-type grid.
	NImportFile::ImportAS3Grid(mConfigContainer.get(), 
			                       mGeometryContainer.get());

	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	std::cout << "Initializing multizone components in: " << std::endl;

	// Loop over each zone and instantiate the necessary objects.
	for(size_t iZone=0; iZone<mSolverContainer.size(); iZone++)
	{
		// Report output.
		std::cout << "  zone: " << iZone << ")" << std::endl;

		// Extract the solver and geometry in this zone.
		auto* zone   = mGeometryContainer->GetZoneGeometry(iZone);
		auto* solver = mSolverContainer[iZone].get();

		// Initialize the physical elements.
		mSolverContainer[iZone]->InitPhysicalElements(mConfigContainer.get(),
				                                          mGeometryContainer.get());
	
		// Initialize the boundary conditions.
		mSolverContainer[iZone]->InitBoundaryConditions(mConfigContainer.get(),
				                                            mGeometryContainer.get());

		// Initialize the solution.
		mInitialContainer->InitializeSolution(mConfigContainer.get(), zone, solver);
	}

	// Report output.
	std::cout << "Done." << std::endl;

	// Display the boundary conditions over all zones.
	NLogger::DisplayBoundaryConditions(mConfigContainer.get(), mSolverContainer);

	// Save initial state, before time-marching.
	mOutputContainer->WriteVisualFile(mConfigContainer.get(), 
			                              mGeometryContainer.get(),
																		mSolverContainer);
}

//-----------------------------------------------------------------------------------

void CDriver::Run
(
 void
)
 /*
	* Function that runs the simulation in time. 
	*/
{
	// Report solver specification as output.
	NLogger::PrintInitSolver(mConfigContainer.get());

	// Extract starting time, ending time and time step.
	const as3double t0 = mConfigContainer->GetStartTime();
	const as3double tf = mConfigContainer->GetFinalTime();
	const as3double dt = mConfigContainer->GetTimeStep();
	
	// Extract the total number of temporal iterations.
	const size_t nIter = mConfigContainer->GetMaxIterTime();


	// Allocate vector for data to monitor.
	as3vector1d<as3double> monitordata(2);


	// Initialize the starting time and iteration count.
	as3double t = t0; size_t i = 0;

	// March in time until either the max iterations or the final time is reached.
	while( (t<tf) && (i<nIter) )
	{

		// Update the solution in time.
		mTemporalContainer->UpdateTime(mConfigContainer.get(),
				                           mGeometryContainer.get(),
																	 mIterationContainer.get(),
																	 mSolverContainer,
																	 t, dt,
																	 monitordata);


		// Update physical time.
		t += dt;

		// Update iteration count.
		i++;

		// Extra processing steps go here.

	}
}




