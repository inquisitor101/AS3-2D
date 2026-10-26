#pragma once 

/*!
 * @file
 * @brief Spatial grid sweeps and element updates within a temporal stage.
 */


#include <cstddef>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "interface_structure.hpp"
#include "openmp_structure.hpp"
#include "collective_residual_computation.hpp"

// Forward declaration to avoid compiler issues.
class CMultizoneSolver;
class ISolver;
class IInterface;
class COpenMP;
class ITemporal;

/*!
 * @brief Spatial residual sweeps and fused element updates for one temporal stage.
 */
class CIteration
{
	public:

		/*!
		 * @brief Constructor of CIteration, which is responsible for an iteration via a single grid
		 *   sweep.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_solver_container solvers and residual storage for all zones.
		 */
		CIteration(const CConfig          *config_container,
				       const CMultizoneSolver *multizone_solver_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CIteration(void);

		/*!
		 * @brief Function that computes the spatial residual and updates the solution for one temporal
		 *   stage.
		 *
		 * @param[in] multizone_geometry_container input geometry container.
		 * @param[in,out] multizone_solver_container input multizone solver container.
		 * @param[in] temporal_container temporal method that updates one element.
		 * @param[in] localtime local physical time.
		 * @param[in] timestep physical time step.
		 * @param[in] iStageRK index of the RK stage.
		 *
		 * @note Called outside a parallel region. This function creates the OpenMP team and
		 *   synchronizes surface writes before element updates.
		 */
		void AdvanceStage(const CMultizoneGeometry *multizone_geometry_container,
		                    CMultizoneSolver         *multizone_solver_container,
		                    const ITemporal          *temporal_container,
		                    as3double                 localtime,
		                    as3double                 timestep,
		                    std::size_t            iStageRK);

		/*!
		 * @brief Function that performs a grid sweep over all the zones.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 * @param[in] openmp_container OpenMP shared parallelization container.
		 * @param[in,out] multizone_solver_container input multizone solver container.
		 * @param[in] localtime local physical time.
		 *
		 * @note Called collectively inside the caller's parallel region. Computes the spatial residual
		 *   without updating the solution.
		 */
		void GridSweep(CConfig            *config_container,
		               CMultizoneGeometry *multizone_geometry_container,
									 COpenMP            *openmp_container,
									 CMultizoneSolver   *multizone_solver_container,
									 as3double           localtime);

	protected:

	private:
		std::size_t mNWorkData = 0; ///< Number of data entries in the work array.

		/*!
		 * @brief Function that computes the residual in all zones.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 * @param[in] openmp_container OpenMP shared parallelization container.
		 * @param[in,out] multizone_solver_container input multizone solver container.
		 * @param[in,out] workarray memory for the working array.
		 * @param[in] localtime local physical time.
		 */
		void ComputeResiduals(CConfig                   *config_container,
		                      CMultizoneGeometry        *multizone_geometry_container,
												  COpenMP                   *openmp_container,
									        CMultizoneSolver          *multizone_solver_container,
									        CPoolMatrixAS3<as3double> &workarray,
												  as3double                  localtime);



		// Disable default constructor.
		CIteration(void) = delete;
		// Disable default copy constructor.
		CIteration(const CIteration&) = delete;
		// Disable default copy operator.
		CIteration& operator=(CIteration&) = delete;
};


