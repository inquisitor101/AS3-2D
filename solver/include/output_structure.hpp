#pragma once 

/*!
 * @file
 * @brief Dispatch of visualization output for the multizone solution.
 */


#include <memory>

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
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 */
		COutput(const CConfig            *config_container,
				    const CMultizoneGeometry *multizone_geometry_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~COutput(void);

		/*!
		 * @brief Function that writes a visualization file.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container geometry of all grid zones.
		 * @param[in] openmp_container OpenMP team configuration.
		 * @param[in] multizone_solver_structure solvers and solution data for all zones.
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
