#pragma once

#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "import_structure.hpp"
#include "output_structure.hpp"
#include "temporal_structure.hpp"
#include "solver_structure.hpp"
#include "initial_condition_structure.hpp"


/*!
 * @brief A class used for the definition, processing and post-processing of all the solver. 
 */
class CDriver
{
	public:
		// Disable default constructor.
		CDriver(void) = delete;
		
		/*!
		 * @brief Constructor of CDriver, which is responsible for the entire solver set-up.
		 */
		CDriver(const char *filename);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CDriver(void);

		/*!
		 * @brief Function that initializes the data for the simulation.
		 */
		void InitializeData(void);

	protected:


	private:
		std::unique_ptr<CConfig>              mConfigContainer;   ///< Container for the configuration options. 
		std::unique_ptr<CGeometry>            mGeometryContainer; ///< Container for the geometry information. 
		std::unique_ptr<COutput>              mOutputContainer;   ///< Container for the output functionalities.
    std::unique_ptr<ITemporal>            mTemporalContainer; ///< Container for the temporal discretization.
		as3vector1d<std::unique_ptr<ISolver>> mSolverContainer;	  ///< Vector of containers for the solver, per each zone.
};
