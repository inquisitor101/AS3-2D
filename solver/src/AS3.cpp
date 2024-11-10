#include "AS3.hpp"

#ifdef HAVE_MPI
	#include "mpi_structure.hpp"
#endif

int main(int argc, char **argv)
{

	// Check whether this is an MPI implementation or not.
#ifdef HAVE_MPI

	// This is an MPI implementation.
	CMPI::Init(&argc, &argv);

	// Output is done on the master rank: 0, only.
	if( CMPI::GetRankI() == 0 )
	{
		NLogger::PrintUsageInstruction(argc, argv);
	}

#else 

	// This is a serial implementation.
	NLogger::PrintUsageInstruction(argc, argv);

#endif

	// Instantiate a driver class and provide it with the configuration file.
	std::unique_ptr<CDriver> driver_container = std::make_unique<CDriver>( argv[1] );

	// Run the simulation.
	driver_container->StartSolver();

	// Exit happily.
	return 0;
}

