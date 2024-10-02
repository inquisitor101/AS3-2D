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
 CConfig                                  *config_container,
 CGeometry                                *geometry_container,
 as3vector1d<std::unique_ptr<IInterface>> &interface_container
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
	for( auto& zone: geometry_container->GetZoneGeometry() )
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

				// Loop over each interface and search for our marker.
				for( auto& interface: interface_container )
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



