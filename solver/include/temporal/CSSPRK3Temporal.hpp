#pragma once

/*!
 * @file
 * @brief Third-order strong-stability-preserving Runge-Kutta time integration.
 */


#include "temporal/ITemporal.hpp"


/*!
 * @brief A class for a temporal discretization based on the SSP-RK3.
 */
class CSSPRK3Temporal final : public ITemporal
{
	public:
	

		/*!
		 * @brief Constructor of CSSPRK3Temporal, which initializes an SSP-RK3 temporal discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CSSPRK3Temporal(const CConfig *config_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CSSPRK3Temporal(void) override;


		/*!
		 * @brief Function that computes the upcoming solution in time, based on an SSP-RK3.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 * @param[in] iteration_container input iteration container.
		 * @param[in] openmp_container OpenMP shared parallelization container.
		 * @param[in,out] solver_container solver and solution storage.
		 * @param[in] physicaltime current physical simulation time.
		 * @param[in] timestep physical time step.
		 */
		void UpdateTime(CConfig            *config_container,
		                CMultizoneGeometry *multizone_geometry_container,
										CIteration         *iteration_container,
										COpenMP            *openmp_container,
										CMultizoneSolver   *solver_container,
										as3double           physicaltime, 
										as3double           timestep) override;

		/*!
		 * @brief Function that applies one SSP-RK3 stage to one element.
		 *
		 * @param[in,out] physical_element solution, residual, and temporal history for one element.
		 * @param[in] timestep physical time step.
		 * @param[in] iStageRK index of the RK stage.
		 *
		 * @note Element-local operation with no mesh traversal or synchronization; the caller owns this
		 *   element for the stage.
		 */
		void UpdateElementSolution(CPhysicalElement *physical_element,
		                             as3double         timestep,
		                             std::size_t    iStageRK) const override;

	protected:

	private:
		std::size_t mNStageRK = 3;     ///< Number of RK evaluation stages.
		as3double      mRk3a[3] = {};  ///< SSP-RK3: a-coefficients.
		as3double      mRk3b[3] = {};  ///< SSP-RK3: b-coefficients.
		as3double      mRk3c[3] = {};  ///< SSP-RK3: c-coefficients.



};
