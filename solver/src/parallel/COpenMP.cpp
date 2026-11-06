#include "parallel/COpenMP.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif


//-----------------------------------------------------------------------------------
// COpenMP member functions.
//-----------------------------------------------------------------------------------


COpenMP::COpenMP
(
 CConfig *config_container
)
 /*
	* Constructor that records the maximum OpenMP team size.
	*/
{
#ifdef HAVE_OPENMP
  mNThread = static_cast<std::size_t>(omp_get_max_threads());
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
	* Destructor for COpenMP.
	*/
{

}





