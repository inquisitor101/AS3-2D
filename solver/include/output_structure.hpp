#pragma once 

#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "vtk_structure.hpp"


/*!
 * @brief A class used for writing information to files. 
 */
class COutput
{
	public:
	
		/*!
		 * @brief Constructor of COutput, which is responsible for the entire output routines.
		 */
		COutput(CConfig   *config_container,
				    CGeometry *geometry_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~COutput(void);

		/*!
		 * @brief Function that writes a visualization file.
		 *
		 * @param[in] config_container pointer to the configuration container.
		 * @param[in] geometry_container pointer to the geometry container.
		 */
		void WriteVisualFile(CConfig                               *config_container,
											   CGeometry                             *geometry_container,
												 as3vector1d<std::unique_ptr<ISolver>> &solver_structure);

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
