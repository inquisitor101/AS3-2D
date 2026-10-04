#include "log_structure.hpp"


//-----------------------------------------------------------------------------------
// NLogger namespace functions.
//-----------------------------------------------------------------------------------


void NLogger::PrintInitSolver
(
 CConfig *config_container
)
 /*
	* Function that prints the information about what solver is used in each zone.
	*/
{
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	
	// Report specific information about solver type and buffer layer, per zone.
	std::cout << "Initiating simulation... " << std::endl;

	// Extract total number of zones.
	auto nZone = config_container->GetnZone();
	
	// Define temporary lambda that prints the type of buffer layer in a zone.
	auto LPrintTypeBufferLayer = [&](ETypeBufferLayer buf) -> std::string
	{
		// String containing the relevant message.
		std::string out;
		
		// Check which type of message to issue, depending on the type of buffer layer.
		switch(buf)
		{
			case(ETypeBufferLayer::NONE): { out = " no buffer layer"; break; }
			default: ERROR("Unkown buffer layer.");
		}

		// return the name of the buffer layer, if any.
		return out;
	};


	// Loop over each zone and print the corresponding solver specs.
	for(unsigned short iZone=0; iZone<nZone; iZone++)
	{
		switch( config_container->GetTypeSolver(iZone) )
		{
			case(ETypeSolver::EE):
			{
	      std::cout << "  Beginning EE Solver in iZone(" << iZone << "): "
					        << LPrintTypeBufferLayer(config_container->GetTypeBufferLayer(iZone));
	    	break;
			}

			default: ERROR("Unknown solver type.");
		}

		std::cout << std::endl;
	}

	std::cout << "Done." << std::endl;
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayBoundaryConditions
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container
)
 /*
	* Function that displays the boundary condition information over all zones.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	std::cout << "Reporting boundary conditions specified:" << std::endl;

	// Max length of marker names, used for output format.
	size_t maxlen = 0;
	
	// Determine the maximum length of the interface marker names.
	for( auto& [name, bc]: config_container->GetMarkerTag() )
	{
		maxlen = std::max( maxlen, name.size() );
	}


	// Loop over each grid zone.
	for( auto& zone: multizone_geometry_container->GetSinglezoneGeometry() )
	{
		// Extract the marker information in this zone.
		auto& marker = zone->GetMarker();

		// Display the zone ID and number of boundaries.
		std::cout << "  zone: " << zone->GetZoneID() << ")\n";
		std::cout << "   ... has " 
			        << marker.size() << " markers with BCs:" << std::endl;

		// Counter for the number of boundaries in this zone.
		size_t n = 0;

		// Loop over each marker and deduce its type and information.
		for( auto& m: marker )
		{
			// Distinguish between regular boundaries and interface/periodic BCs.
			if( m->GetTypeBC() == ETypeBC::INTERFACE )
			{
				// Flag whether the marker is found or not.
				bool found = false;

        // Get the interfaces.
        const auto& interface_container = multizone_solver_container->GetMultizoneInterface();

				// Loop over each interface and search for our marker.
				for( const auto& interface: interface_container )
				{
					// Check if this marker is an owner on this face.
					if( interface->GetIName() == m->GetNameMarker() )
					{
						// Extract zone indices.
						auto  izone = interface->GetIZone();
						auto  jzone = interface->GetJZone();

						// Extract marker names.
						auto& iname = interface->GetIName();
						auto& jname = interface->GetJName();

						// Copy the name tags to a fixed-length string.
						std::string tmp = iname;
						std::string inamepad(maxlen, ' ');
						for(size_t ii=0; ii<tmp.size(); ii++) inamepad[ii] = tmp[ii];

						// Report output concerning the interface markers.
						std::cout << "     "    << izone    << "." 
							        << n++                    << ") INTERFACE: "
											<< " iZone: " << izone    << ", " 
											<< " iName: " << inamepad << " ==> " 
											<< " jZone: " << jzone    << ", "
											<< " jName: " << jname    << std::endl; 
						
						// Flag that we found our marker and break from the loop.
						found = true; break;
					}

					// Check if this marker is a matching pair on this face.
					if( interface->GetJName() == m->GetNameMarker() )
					{
						// Extract zone indices.
						auto  izone = interface->GetJZone();
						auto  jzone = interface->GetIZone();

						// Extract marker names.
						auto& iname = interface->GetJName();
						auto& jname = interface->GetIName();

						// Copy the name tags to a fixed-length string.
						std::string tmp = iname;
						std::string inamepad(maxlen, ' ');
						for(size_t ii=0; ii<tmp.size(); ii++) inamepad[ii] = tmp[ii];

						// Report output concerning the interface markers.
						std::cout << "     "    << izone    << "." 
							        << n++                    << ") INTERFACE: "
											<< " iZone: " << izone    << ", " 
											<< " iName: " << inamepad << " ==> " 
											<< " jZone: " << jzone    << ", "
											<< " jName: " << jname    << std::endl; 
						
						// Flag that we found our marker and break from the loop.
						found = true; break;
					}
				}

				// If we couldn't find the marker, issue an error.
				if( !found ) ERROR("Interface marker could not be located.");
			}
			else
			{
				// This is a regular boundary condition.
				ERROR("Not implemented yet.");
			}
		}
	}


	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void NLogger::MonitorOutput
(
 CConfig       *config_container,
 CMonitorData  *monitor_container,
 unsigned long  iter,
 as3double      time,
 as3double      step
)
 /*
	* Function that outputs the header of the information being displayed.
	*/
{
	// Extract the max number of iterations.
	unsigned long nMaxIter = config_container->GetMaxIterTime();
	// Number of output reports for monitoring progress.
	unsigned long nOutput  = std::max(1ul, nMaxIter/100);
	// Compute number of max digits needed for output.
	unsigned long nDigits  = std::to_string(nMaxIter).size();

	// Display header.
	if( iter%(50*nOutput) == 0 )
	{
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
		std::cout << " Iteration "     << "\t"
			        << " Time "          << "\t" 
							<< " Sync time "     << "\t"
							<< " Steps "         << "\t"
							<< " min(dt) "       << "\t"
							<< " max(dt) "       << "\t"
							<< " Max(Mach) "     << "\n";
		std::cout << "**********************************************"
							<< "**********************************************" << std::endl;
	}

  // Extract the maximum Mach number.
  const as3double Mmax  = monitor_container->mMachMax;
	// Extract the number of substeps per sync.
	const size_t    nsub  = monitor_container->mNSyncSubStep;
	// Extract the min and max time steps, per sync step.
	const as3double dtmin = monitor_container->mMinTimeStep;
	const as3double dtmax = monitor_container->mMaxTimeStep;

	// Display progress.
	std::cout << std::scientific 
						<< std::setprecision(6)
		        << " " 
						<< std::setw(static_cast<int>(nDigits)) << iter
						<< "\t"  << time
						<< "\t"  << step
						<< "\t"  << nsub
						<< "\t"  << dtmin
						<< "\t"  << dtmax
            << "\t"  << Mmax
						<< std::endl;
}

//-----------------------------------------------------------------------------------

void NLogger::DisplayOpenMPInfo
(
 COpenMP                  *openmp_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const CMultizoneSolver   *multizone_solver_container
)
 /*
	* Function that displays the OpenMP information, if any.
	*/
{
  // Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
#ifdef HAVE_OPENMP
  // Get max number of threads specified.
  const size_t nThreads = omp_get_max_threads();
  std::cout << "This is a parallel implementation using: "
            << nThreads << " threads." << std::endl;

	// Get the total number of elements in all zones.
	const size_t nElemTotal = multizone_geometry_container->GetnElemTotal();
	// Get the total number of i-faces in all zones.
	const size_t nIFace     = multizone_geometry_container->GetnIFace();
	// Get the total number of j-faces in all zones.
	const size_t nJFace     = multizone_geometry_container->GetnJFace();

  // Estimate computational work load of each thread.
  as3vector1d<size_t> workloadDOFs(nThreads, 0);
	as3vector1d<size_t> workloadElem(nThreads, 0);
	as3vector1d<size_t> workloadIFace(nThreads, 0);
	as3vector1d<size_t> workloadJFace(nThreads, 0);

  as3vector1d<size_t> workloadIFaceInternal_standard(nThreads,  0);
  as3vector1d<size_t> workloadIFaceBoundary_standard(nThreads,  0);
  as3vector1d<size_t> workloadIFaceInterface_standard(nThreads, 0);

  as3vector1d<size_t> workloadIFaceInternal_balanced(nThreads,  0);
  as3vector1d<size_t> workloadIFaceBoundary_balanced(nThreads,  0);
  as3vector1d<size_t> workloadIFaceInterface_balanced(nThreads, 0);

  as3vector1d<size_t> workloadJFaceInternal_standard(nThreads,  0);
  as3vector1d<size_t> workloadJFaceBoundary_standard(nThreads,  0);
  as3vector1d<size_t> workloadJFaceInterface_standard(nThreads, 0);

  as3vector1d<size_t> workloadJFaceInternal_balanced(nThreads,  0);
  as3vector1d<size_t> workloadJFaceBoundary_balanced(nThreads,  0);
  as3vector1d<size_t> workloadJFaceInterface_balanced(nThreads, 0);

	// Estimate the i-surface workload.
#pragma omp parallel for schedule(static)
	for(size_t i=0; i<nIFace; i++)
	{
    // Thread index.
    const size_t iThread = omp_get_thread_num();

		// Accumulate the number of elements per thread.
		workloadIFace[iThread]++;

    // Get the relevant face information for the standard and load-balanced partitions.
    const auto face_info_standard = multizone_geometry_container->GetFlattenedIndexIFace(i);
    const auto face_info_balanced = multizone_geometry_container->GetFlattenedIndexIFaceLoadBalanced(i);

    // Estimate the work load for the standard approach.
    switch( face_info_standard.mFaceType )
    {
      case( ETypeFaceGeometry::INTERNAL ):  { workloadIFaceInternal_standard[iThread]++;  break; }
      case( ETypeFaceGeometry::BOUNDARY ):  { workloadIFaceBoundary_standard[iThread]++;  break; }
      case( ETypeFaceGeometry::INTERFACE ): { workloadIFaceInterface_standard[iThread]++; break; }
      default: ERROR("Unknown face type.");
    }

    // Then, estimate the workload for the load-balanced approach.
    switch( face_info_balanced.mFaceType )
    {
      case( ETypeFaceGeometry::INTERNAL ):  { workloadIFaceInternal_balanced[iThread]++;  break; }
      case( ETypeFaceGeometry::BOUNDARY ):  { workloadIFaceBoundary_balanced[iThread]++;  break; }
      case( ETypeFaceGeometry::INTERFACE ): { workloadIFaceInterface_balanced[iThread]++; break; }
      default: ERROR("Unknown face type.");
    }
	}

	// Estimate the j-surface workload.
#pragma omp parallel for schedule(static)
	for(size_t i=0; i<nJFace; i++)
	{
    // Thread index.
    const size_t iThread = omp_get_thread_num();

		// Accumulate the number of elements per thread.
		workloadJFace[iThread]++;
	
    // Get the relevant face information for the standard and load-balanced partitions.
    const auto face_info_standard = multizone_geometry_container->GetFlattenedIndexJFace(i);
    const auto face_info_balanced = multizone_geometry_container->GetFlattenedIndexJFaceLoadBalanced(i);

    // Estimate the work load for the standard approach.
    switch( face_info_standard.mFaceType )
    {
      case( ETypeFaceGeometry::INTERNAL ):  { workloadJFaceInternal_standard[iThread]++;  break; }
      case( ETypeFaceGeometry::BOUNDARY ):  { workloadJFaceBoundary_standard[iThread]++;  break; }
      case( ETypeFaceGeometry::INTERFACE ): { workloadJFaceInterface_standard[iThread]++; break; }
      default: ERROR("Unknown face type.");
    }

    // Then, estimate the workload for the load-balanced approach.
    switch( face_info_balanced.mFaceType )
    {
      case( ETypeFaceGeometry::INTERNAL ):  { workloadJFaceInternal_balanced[iThread]++;  break; }
      case( ETypeFaceGeometry::BOUNDARY ):  { workloadJFaceBoundary_balanced[iThread]++;  break; }
      case( ETypeFaceGeometry::INTERFACE ): { workloadJFaceInterface_balanced[iThread]++; break; }
      default: ERROR("Unknown face type.");
    }
  }

	// Estimate the elements' workload.
#pragma omp parallel for schedule(static)
  for(size_t i=0; i<nElemTotal; i++)
	{
		// Extract the element indices.
		const auto elem_info = multizone_geometry_container->GetFlattenedIndexVolumeElement(i);

		// Deduce the current element's zone and index.
		const auto iZone = elem_info.mIndexZone; 
		const auto iElem = elem_info.mIndexElem; 

    // Thread index.
    const size_t iThread = omp_get_thread_num();

		// Accumulate the number of elements per thread.
		workloadElem[iThread]++;

    // Accumulate the number of solution DOFs per thread.
		workloadDOFs[iThread] += multizone_solver_container->GetSinglezoneSolver(iZone)->GetStandardElement()->GetnSol2D(); 
  }


  // Compute the number of max digits needed for the output.
  size_t nDigits = 0;
  for( auto& work: workloadDOFs ) nDigits = std::max( nDigits, work );
  // Deduce the max value needed for the digits width.
  nDigits = std::to_string(nDigits).size();

  // Total number of work load.
  std::cout << "The estimated workload shared among each thread is:\n";
  for(size_t i=0; i<workloadDOFs.size(); i++)
    std::cout << "  Thread(" << i << ") has:\n" 
			        << "   (*) " << std::setw(nDigits)
			        << workloadIFace[i] << " [nIFace/thread]\n"
							<< "      -> " << std::setw(nDigits)
              << workloadIFaceInternal_standard[i]  << " [nInternalIFace/thread]  (standard)\t" << std::setw(nDigits)  
              << workloadIFaceInternal_balanced[i]  << " [nInternalIFace/thread]  (balanced)\n" << std::setw(nDigits)
              << "      -> " << std::setw(nDigits)
              << workloadIFaceBoundary_standard[i]  << " [nBoundaryIFace/thread]  (standard)\t" << std::setw(nDigits)
              << workloadIFaceBoundary_balanced[i]  << " [nBoundaryIFace/thread]  (balanced)\n" << std::setw(nDigits)
              << "      -> " << std::setw(nDigits)
              << workloadIFaceInterface_standard[i] << " [nInterfaceIFace/thread] (standard)\t" << std::setw(nDigits)
              << workloadIFaceInterface_balanced[i] << " [nInterfaceIFace/thread] (balanced)\n" << std::setw(nDigits)
              << "   (*) " << std::setw(nDigits)
							<< workloadJFace[i] << " [nJFace/thread]\n"
							<< "      -> " << std::setw(nDigits)
              << workloadJFaceInternal_standard[i]  << " [nInternalJFace/thread]  (standard)\t" << std::setw(nDigits) 
              << workloadJFaceInternal_balanced[i]  << " [nInternalJFace/thread]  (balanced)\n" << std::setw(nDigits)
              << "      -> " << std::setw(nDigits)
              << workloadJFaceBoundary_standard[i]  << " [nBoundaryJFace/thread]  (standard)\t" << std::setw(nDigits)
              << workloadJFaceBoundary_balanced[i]  << " [nBoundaryJFace/thread]  (balanced)\n" << std::setw(nDigits)
              << "      -> " << std::setw(nDigits)
              << workloadJFaceInterface_standard[i] << " [nInterfaceJFace/thread] (standard)\t" << std::setw(nDigits)
              << workloadJFaceInterface_balanced[i] << " [nInterfaceJFace/thread] (balanced)\n" << std::setw(nDigits)
							<< "   (*) " << std::setw(nDigits)
			        << workloadElem[i] << " [nElement/thread]\n"
							<< "   (*) " << std::setw(nDigits)
              << workloadDOFs[i] << " [nSolDOFs/thread]" << std::endl;
#else
	std::cout << "This is a serial implementation." << std::endl;
#endif
}








