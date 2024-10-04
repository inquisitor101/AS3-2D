#include "output_structure.hpp"


//-----------------------------------------------------------------------------------
// COutput member functions.
//-----------------------------------------------------------------------------------

COutput::COutput
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Constructor for the output class, which is responsible for the entire output routines.
	*/
{
	mVTKContainer = CGenericFactory::CreateVTKContainer(config_container, geometry_container);
}

//-----------------------------------------------------------------------------------

COutput::~COutput
(
 void
)
 /*
	* Destructor, which cleans up after the output class.
	*/
{

}

//-----------------------------------------------------------------------------------

void COutput::WriteVisualFile
(
 CConfig                               *config_container,
 CGeometry                             *geometry_container,
 COpenMP                               *openmp_container,
 as3vector1d<std::unique_ptr<ISolver>> &solver_container
)
 /*
	* Function that writes a visualization file.
	*/
{
	// Ensure the VTK container is initialized.
	if( mVTKContainer ) mVTKContainer->WriteFileVTK(config_container, 
			                                            geometry_container,
																									openmp_container,
																									solver_container);
}
