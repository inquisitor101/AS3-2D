#include "openmp_structure.hpp"


//-----------------------------------------------------------------------------------
// COpenMP member functions.
//-----------------------------------------------------------------------------------


COpenMP::COpenMP
(
 CConfig *config_container
)
 /*
	* Constructor for the OpenMP shared memory parallelization class.
	*/
{
#ifdef HAVE_OPENMP
  mNThread = omp_get_max_threads();
#else
  mNThread = 1;
#endif
}

//-----------------------------------------------------------------------------------

COpenMP::~COpenMP
(
 void
)
 /*
	* Destructor, which cleans up after the OpenMP class.
	*/
{

}





