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

	// Initialize the data monitoring container.
	mMonitoringContainer = CGenericFactory::CreateMonitoringContainer(mConfigContainer.get());
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

	// Extract the interface boundaries.
	auto& interface = mConfigContainer->GetInterfaceParamMarker();

	// If interface conditions are specified, initialize them.
	if( interface.size() )
	{
		// Allocate the required number of interface containers.
		mInterfaceContainer.reserve( interface.size() );

		// Initialize the interface boundaries.
		for(size_t i=0; i<interface.size(); i++)
		{
			mInterfaceContainer.emplace_back
			(
			 CGenericFactory::CreateInterfaceContainer(mConfigContainer.get(), 
					                                       mGeometryContainer.get(),
																								 interface[i].get(),
																								 mSolverContainer)
			);
		}
	}

	// Initialize the OpenMP container.
	mOpenMPContainer = std::make_unique<COpenMP>(mConfigContainer.get(), 
			                                         mGeometryContainer.get(),
																							 mSolverContainer);

	// Report output.
	std::cout << "Done." << std::endl;

	// Display the boundary conditions over all zones.
	NLogger::DisplayBoundaryConditions(mConfigContainer.get(), 
			                               mGeometryContainer.get(), 
																		 mInterfaceContainer);
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
	// For convenience, extract the properties of the simulation end time.
	const size_t nIter = mConfigContainer->GetMaxIterTime();
	// Extract the visualization output frequency.
	const size_t fvis  = mConfigContainer->GetWriteVisFreq();

	// Flag whether the the visualization file is written.
	bool isvis = false;

	// Check if we need to write the visualization file.
	if( i%fvis == 0 )
	{
		mOutputContainer->WriteVisualFile(mConfigContainer.get(), 
				                              mGeometryContainer.get(),
																			mOpenMPContainer.get(),
																			mSolverContainer);
	
		// Update the visualization flag.
		isvis = true;
	}

	// Check if this is the end of the simulation.
	if( i >= nIter )
	{

		// Write the visualization file, if it hasnt been written.
		if( !isvis )
		{
			mOutputContainer->WriteVisualFile(mConfigContainer.get(), 
					                              mGeometryContainer.get(),
																				mOpenMPContainer.get(),
																				mSolverContainer);
		}
	}

	// Monitor data.
	NLogger::MonitorOutput(mConfigContainer.get(), 
			                   mMonitoringContainer.get(),
				                 i, t, dt);

	// TODO: write GNUplot file

  // Check for floating-point errors at run-time.
#ifdef ENABLE_NAN_CHECK
	NError::CheckFloatingError();
#endif
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

	// Extract the simulation's starting time.
	const as3double t0 = mConfigContainer->GetStartTime();

	// Extract the synchronization time step.
	const as3double ts = mConfigContainer->GetTimeStep();

	// Extract the total number of time synchronization iterations.
	const size_t nSyncStep = mConfigContainer->GetMaxIterTime();

	// Initialize the current time and iteration step.
	as3double cTime = t0; size_t iSyncStep = 0;

	// Write the initial output.
	WriteOutput(iSyncStep, cTime, ts);


	// March in time until the max number of iterations is reached.
	while( iSyncStep < nSyncStep )
	{
		// Compute the solution for the current synchronized step.
		ExecuteTimeSyncStep(cTime);

		// Update physical time and iteration count.
		cTime += ts; iSyncStep++;
	
		// Write the output data, if need be.
		WriteOutput(iSyncStep, cTime, ts);
	}
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
	// Extract the specified CFL number.
	const as3double cfl = mConfigContainer->GetCFL();

	// Initialize the inverse of the lowest stable time step.
	as3double dtinv = C_ZERO;

	// Max Mach number squared, used for monitoring.
	as3double maxM2 = C_ZERO;

	// Get the total number of elements in all zones.
	const size_t nElemTotal = mOpenMPContainer->GetnElemTotal();

	// Loop over all the elements in all the solvers.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), reduction(max:dtinv), reduction(max:maxM2)
#endif
	for(size_t i=0; i<nElemTotal; i++)
	{
		// Deduce the current element's zone and index.
		const unsigned short iZone = mOpenMPContainer->GetIndexVolume(i)->mZone;
		const unsigned int   iElem = mOpenMPContainer->GetIndexVolume(i)->mElem;

		// Extract the relevant solver.
		auto& solver  = mSolverContainer[iZone];
		// Extract the relevant physical element.
		auto* element = solver->GetPhysicalElement(iElem);

		// Extract the number of DOFs on each element in this zone. 
		const size_t   nSol2D = solver->GetStandardElement()->GetnSol2D();
		// Extract the polynomial order in this zone.
		const as3double npoly = static_cast<as3double>( solver->GetStandardElement()->GetnPolySol() );
		// Deduce the maximum inviscid polynomial coefficient.
		const as3double f1    = npoly*npoly;

		// Extract the solution.
		const auto& sol = element->mSol2D;	

		// Extract average normals, based on the surface directions (i and j).
		const as3double *ni = element->mAvgNormIDir;
		const as3double *nj = element->mAvgNormJDir;

		// Compute the inverse of the average length scale in the i and j-direction.
		const as3double ovli = C_ONE/element->mLengthScaleIDir;
		const as3double ovlj = C_ONE/element->mLengthScaleJDir;

		// Initialize the max of the inverse (inviscid) spectral radius on this element.
		as3double maxsrinv = C_ZERO;
		// Initialize the max of the Mach squared on this element. 
		as3double maxm2    = C_ZERO;

		// Loop over each DOF and compute the stability time limit.
		for(size_t l=0; l<nSol2D; l++)
		{
  		// Compute the primitive variables.
  		const as3double rho   = sol(0,l);
  		const as3double ovrho = C_ONE/rho;
  		const as3double u     = ovrho*sol(1,l);
  		const as3double v     = ovrho*sol(2,l);
  		const as3double p     = C_GM1*( sol(3,l) - C_HALF*(u*sol(1,l) + v*sol(2,l)) );

			// Compute the speed of sound  and its squared.
			const as3double a2    = C_GMA*p*ovrho;
			const as3double a     = std::sqrt(a2);

			// Compute the average i and j-projected velocities.
			const as3double ui    = u*ni[0] + v*ni[1];
			const as3double uj    = u*nj[0] + v*nj[1];

			// Compute the maximum eigenvalues of the inviscid terms (which are acoustic).
			const as3double lmbi  = std::abs(ui) + a;
			const as3double lmbj  = std::abs(uj) + a;

			// Compute the inverse of the spectral radius of the inviscid terms.
			const as3double srinv = lmbi*ovli + lmbj*ovlj; 
			
			// Compute the max of the spectral radius inverted.
			maxsrinv  = std::max( srinv, maxsrinv );
		
    	// Compute the local Mach number squared.
    	const as3double M2 = (u*u + v*v)/a2;
    	// Check if this value is the largest.
    	maxm2 = std::max(maxm2, M2);

			// Ensure the speed of sound is positive.
			if( a2 < C_ZERO ) ERROR("Negative speed of sound encountered.");
		}

		// Assign the max of the inverse time step.
		dtinv = std::max( f1*maxsrinv, dtinv );
		// Assign the max of the Mach squared.
		maxM2 = std::max( maxm2, maxM2 );
	}

	// Update the monitored data.
	mMonitoringContainer->mMachMax = std::sqrt(maxM2);

	// Estimate the stable inverse time step.
	const as3double dt = cfl/dtinv;

	// Return the expected total number of synchronization time steps.
	return dt;
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
	// Flag that specified whether the time is sync'd or not.
	bool issync = false;

	// Counter for the number of sub-steps in this function.
	size_t nSubStep = 0;

	// Minimum and maximum time steps.
	as3double dtmin = static_cast<as3double>( 99999.0 );
	as3double dtmax = C_ZERO;

	// Extract synchronization time step.
	const as3double ts = mConfigContainer->GetTimeStep();
	// Deduce the final time, in this synchronization time step.
	const as3double tf = t0 + ts;
	// Cutoff time, which adjusts the remaining time step to reach synchronization.
	const as3double tc = t0 + static_cast<as3double>(0.99)*ts;

	// Current elapsed time.
	as3double time = t0;

	// Loop until the physical time is synchronized.
	while( !issync )
	{
		// Compute the time step.
		as3double dt = ComputeTimeStep();

		// If the stable time step is larger than the synchronization step, issue a warning.
		if( dt > ts )
		{
			WARNING("Inefficient synchronization time step: dt(sync)/dt(stable) = " + std::to_string((dt/ts)));

			// Set the stable time step to the sync step.
			dt = std::min(dt, ts);
		}
		
		// Check whether the synchronization time is reached.
		if( (time+dt) > tc ) 
		{
			// Set the remaining time step and flag that sync is reached.
			dt = tf - time; issync = true;
		}


		// Update the solution in time.
		mTemporalContainer->UpdateTime(mConfigContainer.get(),
				                           mGeometryContainer.get(),
																	 mIterationContainer.get(),
																	 mOpenMPContainer.get(),
																	 mSolverContainer,
																	 mInterfaceContainer,
																	 time, dt);

		// Update the elapsed time and increment the number of sub-steps.
		time += dt; nSubStep++;

		// Compute the global min and max time steps, per synchronization.
		dtmin = std::min( dt, dtmin );
		dtmax = std::max( dt, dtmax );
	}

	// Book-keep the number of substeps, min and max time steps.
	mMonitoringContainer->mNSyncSubStep = nSubStep;
	mMonitoringContainer->mMinTimeStep  = dtmin;
	mMonitoringContainer->mMaxTimeStep  = dtmax;
}





