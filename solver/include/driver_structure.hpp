#pragma once

/*!
 * @file
 * @brief Simulation setup, synchronized time advancement and output.
 */


#include <cstddef>
#include <memory>

#include <ctime>
#include <chrono>

#include "log_structure.hpp"
#include "option_structure.hpp"
#include "config_structure.hpp"
#include "import_structure.hpp"
#include "output_structure.hpp"
#include "solver_structure.hpp"
#include "openmp_structure.hpp"
#include "factory_structure.hpp"
#include "geometry_structure.hpp"
#include "temporal_structure.hpp"
#include "iteration_structure.hpp"
#include "interface_structure.hpp"
#include "monitoring_structure.hpp"
#include "initial_condition_structure.hpp"


/*!
 * @brief Simulation construction, initialization, time advancement and output.
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
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CDriver(void);

		/*!
		 * @brief Function that is responsible for the entire simulation, from start to end.
		 */
		void StartSolver(void);

	protected:

	private:
		std::unique_ptr<CConfig>                 mConfigContainer;            ///< Container for the configuration options. 	
		std::unique_ptr<COutput>                 mOutputContainer;            ///< Container for the output functionalities.
    std::unique_ptr<CMonitorData>            mMonitoringContainer;        ///< Container for the monitoring data.
		std::unique_ptr<COpenMP>                 mOpenMPContainer;            ///< Container for the OpenMP functionalities.
		std::unique_ptr<ITemporal>               mTemporalContainer;          ///< Container for the temporal discretization.
		std::unique_ptr<CIteration>              mIterationContainer;         ///< Container for a single grid-sweep iteration.
		std::unique_ptr<IInitialCondition>       mInitialContainer;           ///< Container for the initial condition.
    std::unique_ptr<CMultizoneGeometry>      mMultizoneGeometryContainer; ///< Container for the multizone geometries. 
    std::unique_ptr<CMultizoneSolver>        mMultizoneSolverContainer;   ///< Container for the multizone solvers.

		/*!
		 * @brief Function that advances the simulation through its configured synchronization
		 *   intervals.
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

		/*!
		 * @brief Function that writes the output information.
		 *
		 * @param[in] i temporal iteration count.
		 * @param[in] t current physical time.
		 * @param[in] dt current time step.
		 */
		void WriteOutput(std::size_t i, 
				             as3double   t, 
										 as3double   dt);

		/*!
		 * @brief Function that computes the stable global time step from the CFL constraint and element
		 *   states.
		 *
		 * @return stable time step.
		 */
		as3double ComputeTimeStep(void);

		/*!
		 * @brief Function that computes the solution based on a synchronized time step.
		 *
		 * @param[in] t0 current starting time in the synchronization.
		 */
		void ExecuteTimeSyncStep(as3double t0);

		// Disable default constructor.
		CDriver(void) = delete;
		// Disable default copy constructor.
		CDriver(const CDriver&) = delete;
		// Disable default copy operator.
		CDriver& operator=(CDriver&) = delete;
};
