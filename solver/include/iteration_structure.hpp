#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"


// Forward declaration to avoid compiler issues.
class ISolver;


/*!
 * @brief A class for a single grid sweep across all zones. 
 */
class CIteration
{
	public:

		/*!
		 * @brief Constructor of CIteration, which is responsible for an iteration via a single grid sweep.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CIteration(CConfig                               *config_container,
				       as3vector1d<std::unique_ptr<ISolver>> &solver_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CIteration(void);

		/*!
		 * @brief Function that performs a grid sweep over all the zones.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] localtime local physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void GridSweep(CConfig                               *config_container,
		               CGeometry                             *geometry_container,
									 as3vector1d<std::unique_ptr<ISolver>> &solver_container,
									 as3double                              localtime, 
									 as3vector1d<as3double>                &monitordata);

	protected:

	private:
		size_t mNWorkData = 0; // Number of data entries in the work array.

		/*!
		 * @brief Function that preprocesses the solution, before sweeping the grid.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] localtime local physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void PreProcessIteration(CConfig                               *config_container,
		                         CGeometry                             *geometry_container,
									           as3vector1d<std::unique_ptr<ISolver>> &solver_container,
									           as3double                              localtime, 
														 CPoolMatrixAS3<as3double>             &workarray,
									           as3vector1d<as3double>                &monitordata);

		/*!
		 * @brief Function that postprocesses the solution, after sweeping the grid.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] localtime local physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void PostProcessIteration(CConfig                               *config_container,
		                          CGeometry                             *geometry_container,
									            as3vector1d<std::unique_ptr<ISolver>> &solver_container,
									            as3double                              localtime, 
															CPoolMatrixAS3<as3double>             &workarray,
									            as3vector1d<as3double>                &monitordata);

		/*!
		 * @brief Function that computes the residual in all zones.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] localtime local physical time.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void ComputeResidual(CConfig                               *config_container,
		                     CGeometry                             *geometry_container,
									       as3vector1d<std::unique_ptr<ISolver>> &solver_container,
									       as3double                              localtime,
												 CPoolMatrixAS3<as3double>             &workarray,
									       as3vector1d<as3double>                &monitordata);



		// Disable default constructor.
		CIteration(void) = delete;
		// Disable default copy constructor.
		CIteration(const CIteration&) = delete;
		// Disable default copy operator.
		CIteration& operator=(CIteration&) = delete;
};


