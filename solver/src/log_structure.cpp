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
	for(auto iZone=0; iZone<nZone; iZone++)
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
 CConfig                               *config_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
 /*
	* Function that displays the boundary condition information over all zones.
	*/
{
	// Report output.
	std::cout << "----------------------------------------------"
							 "----------------------------------------------\n";
	std::cout << "Reporting boundary conditions specified:" << std::endl;


	// Loop over all the zones.
	for( auto& solver: solver_container )
	{
		// Extract the boundary container in this zone.
		auto& boundary_container = solver->GetBoundaryContainer();

		// Display the zone ID and number of boundaries.
		std::cout << "  zone: " << solver->GetZoneID() << ")\n";
		std::cout << "   ... has " 
			        << boundary_container.size() << " markers with BCs:" << std::endl;

		// Loop over all the boundaries and print their types.
		size_t n = 0;
		for( auto& bc : boundary_container )
		{
			// If this is a periodic boundary, treat it differently.
			if( bc->GetMarker()->GetMarkerBC() == ETypeBC::PERIODIC )
			{

				// Max length of marker names, used for output format.
				size_t maxlen = 0;
	
				// Determine the maximum length of the periodic marker names.
				for( auto& name: config_container->GetMarkerNamePeriodic() )
				{
					maxlen = std::max( maxlen, name[0].size() );
				}

				// Cast the boundary to a periodic boundary.
				CPeriodicBoundary *pbc = dynamic_cast<CPeriodicBoundary*>( bc.get() );
				// This must be succesful, otherwise, Housten, we have a problem.
				if( !pbc ) ERROR("Boundary condition must be periodic.");

				// Abbreviate the pointers to the current and matching markers.
				auto* im = bc->GetMarker();
				auto* jm = pbc->GetMatchingMarker(); 

				// Copy the name tags to a fixed-length string.
				std::string tmp = im->GetNameMarkerTag();
				std::string iMarkPadded(maxlen, ' ');
				for(size_t ii=0; ii<tmp.size(); ii++) iMarkPadded[ii] = tmp[ii];

				// Report a meaningful output.
				std::cout << "     " << solver->GetZoneID()       << "." 
					        << n << ") PERIODIC: "
									<< " iZone: " << im->GetZoneID()        << ", " 
									<< " iName: " << iMarkPadded            << " ==> " 
									<< " jZone: " << jm->GetZoneID()        << ", "
									<< " jName: " << jm->GetNameMarkerTag() << std::endl; 
			}
			else
			{
				// For now, issue an error, as nothing is implemented.
				ERROR("non-periodic BCs have not been implemented (yet).");
			}

			// Increment the number of boundaries.
			n++;
		}
	}


	// Report output.
	std::cout << "Done." << std::endl;
}





