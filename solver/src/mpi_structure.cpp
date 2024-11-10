#include "mpi_structure.hpp"


//-----------------------------------------------------------------------------------
// CMPI member functions.
//-----------------------------------------------------------------------------------


CMPI::CMPI
(
 CConfig *config_container
)
 /*
	* Constructor for the MPI (distributed memory) parallelization class.
	*/
{
	// Get the current rank and the total number of ranks specified.
	// First of all, ensure that the number of zones matches the number of ranks.
	if( config_container->GetnZone() != mNRank )
	{
		ERROR("Number of ranks must be equal to the number of zones, i.e. 1 zone/rank.");
	}
}


CMPI::~CMPI
(
 void
)
 /*
	* Destructor, which cleans up after the MPI class.
	*/
{
	MPI_Finalize();
}


