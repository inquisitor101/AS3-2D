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
}

