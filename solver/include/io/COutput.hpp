#pragma once

/*!
 * @file
 * @brief Dispatch of visualization output for the multizone solution.
 */


#include <memory>

#include "common/AS3ScalarTypes.hpp"

class CConfig;
class CMultizoneGeometry;
class CMultizoneSolver;
class COpenMP;
class IFileVTK;


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
		 * @param[in] time physical simulation time of the snapshot [s].
		 */
		void WriteVisualFile(const CConfig            *config_container,
		                     const CMultizoneGeometry *multizone_geometry_container,
		                     const COpenMP            *openmp_container,
		                     const CMultizoneSolver   *multizone_solver_structure,
		                     as3double                 time);

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
