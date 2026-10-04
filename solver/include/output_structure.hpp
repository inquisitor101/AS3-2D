#pragma once 

#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "vtk_structure.hpp"
#include "openmp_structure.hpp"


/*!
 * @brief A class used for writing information to files. 
 */
class COutput
{
	public:
	
		/*!
		 * @brief Constructor of COutput, which is responsible for the entire output routines.
		 */
		COutput(const CConfig            *config_container,
				    const CMultizoneGeometry *multizone_geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~COutput(void);

		/*!
		 * @brief Function that writes a visualization file.
		 */
		void WriteVisualFile(const CConfig            *config_container,
											   const CMultizoneGeometry *multizone_geometry_container,
												 const COpenMP            *openmp_container,
												 const CMultizoneSolver   *multizone_solver_structure);

	protected:

	private:
		std::unique_ptr<IFileVTK> mVTKContainer;  ///< Container for a VTK file format.

		// Disable default constructor.
		COutput(void) = delete;
		// Disable default copy constructor.
		COutput(const COutput&) = delete;
		// Disable default copy operator.
		COutput& operator=(COutput&) = delete;	
};
