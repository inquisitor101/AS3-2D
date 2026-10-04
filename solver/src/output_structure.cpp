#include "output_structure.hpp"


//-----------------------------------------------------------------------------------
// COutput member functions.
//-----------------------------------------------------------------------------------

COutput::COutput
(
 const CConfig   *config_container,
 const CGeometry *geometry_container
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
 const CConfig          *config_container,
 const CGeometry        *geometry_container,
 const COpenMP          *openmp_container,
 const CMultizoneSolver *multizone_solver_container
)
 /*
	* Function that writes a visualization file.
	*/
{
	// Ensure the VTK container is initialized.
	if( mVTKContainer ) mVTKContainer->WriteFileVTK(config_container, 
			                                            geometry_container,
																									openmp_container,
																									multizone_solver_container);
}
