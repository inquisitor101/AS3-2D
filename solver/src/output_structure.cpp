#include "output_structure.hpp"


//-----------------------------------------------------------------------------------
// COutput member functions.
//-----------------------------------------------------------------------------------

COutput::COutput
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Constructor that creates the visualization writer.
	*/
{
	mVTKContainer = CGenericFactory::CreateVTKContainer(config_container, multizone_geometry_container);
}

//-----------------------------------------------------------------------------------

COutput::~COutput
(
 void
)
 /*
	* Destructor for COutput.
	*/
{

}

//-----------------------------------------------------------------------------------

void COutput::WriteVisualFile
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container,
 const COpenMP            *openmp_container,
 const CMultizoneSolver   *multizone_solver_container
)
 /*
	* Function that writes the current multizone solution through the visualization writer.
	*/
{
	// Ensure the VTK container is initialized.
	if( mVTKContainer ) mVTKContainer->WriteFileVTK(config_container, 
			                                            multizone_geometry_container,
																									openmp_container,
																									multizone_solver_container);
}
