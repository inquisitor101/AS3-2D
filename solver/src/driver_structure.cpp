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
	const std::chrono::duration<as3double> lapsedtime_phys = phys_t1 - phys_t0;

	// Report lapsed time.
	std::cout << "\n% % % % % % % % % % % % % % % % % % % % % % %" << std::endl;
	std::cout << std::scientific << std::setprecision(10) << std::setw(10)
	          << "lapsed (proc) time [sec]: " << lapsedtime_proc << "s %" << "\n"
	          << "lapsed (phys) time [sec]: " << lapsedtime_phys <<  " %" << std::endl;
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

}

//-----------------------------------------------------------------------------------

void CDriver::WriteOutput
(
 size_t    i,
 as3double t,
 as3double dt
)
 /*
	* Function that writes the output information. 
	*/
{

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

}

//-----------------------------------------------------------------------------------

as3double CDriver::ComputeTimeStep
(
 void
)
 /*
	* Function that computes the time step, based on stability constraints. 
	*/
{

	return 0.0;
}

//-----------------------------------------------------------------------------------

void CDriver::ExecuteTimeSyncStep
(
 as3double t0
)
 /*
	* Function that computes the solution based on a synchronized time step.
	*/
{

}





