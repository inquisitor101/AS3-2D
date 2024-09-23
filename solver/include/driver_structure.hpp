#pragma once

#include <ctime>
#include <chrono>

#include "log_structure.hpp"
#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "import_structure.hpp"
#include "output_structure.hpp"
#include "temporal_structure.hpp"
#include "solver_structure.hpp"
#include "iteration_structure.hpp"
#include "interface_structure.hpp"
#include "initial_condition_structure.hpp"


/*!
 * @brief A class used for the definition, processing and post-processing of all the solver. 
 */
class CDriver
{
	public:

		/*!
		 * @brief Constructor of CDriver, which is responsible for the entire solver set-up.
		 *
		 * @param[in] filename input configuration filename.
		 */
		CDriver(const char *filename);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CDriver(void);

		/*!
		 * @brief Function that is responsible for the entire simulation, from start to end.
		 */
		void StartSolver(void);

	protected:

	private:
		std::unique_ptr<CConfig>                 mConfigContainer;    ///< Container for the configuration options. 
		std::unique_ptr<CGeometry>               mGeometryContainer;  ///< Container for the geometry information. 
		std::unique_ptr<COutput>                 mOutputContainer;    ///< Container for the output functionalities.
    std::unique_ptr<ITemporal>               mTemporalContainer;  ///< Container for the temporal discretization.
		std::unique_ptr<CIteration>              mIterationContainer; ///< Container for a single grid-sweep iteration.
		std::unique_ptr<IInitialCondition>       mInitialContainer;   ///< Container for the initial condition.
		as3vector1d<std::unique_ptr<ISolver>>    mSolverContainer;	  ///< Vector of containers for the solver, per each zone.
		as3vector1d<std::unique_ptr<IInterface>> mInterfaceContainer; ///< Vector of containers for the zone interfaces.

		/*!
		 * @brief Function that marches the solution in time.
		 */
		void Run(void);

		/*!
		 * @brief Function that preprocesses the data, prior to starting the simulation.
		 */
		void PreProcess(void);

		/*!
		 * @brief Function that initializes the data for the simulation.
		 */
		void InitializeData(void);


		// Disable default constructor.
		CDriver(void) = delete;
		// Disable default copy constructor.
		CDriver(const CDriver&) = delete;
		// Disable default copy operator.
		CDriver& operator=(CDriver&) = delete;

};
