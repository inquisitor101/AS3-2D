#pragma once

/*!
 * @file
 * @brief Time integration interface.
 */


#include <cstddef>

#include "common/AS3Options.hpp"

class CConfig;
class CMultizoneGeometry;
class CMultizoneSolver;
class CIteration;
class COpenMP;
class CPhysicalElement;


/*!
 * @brief An interface class used for the temporal discretization.
 */
class ITemporal
{
	public:

		/*!
		 * @brief Constructor of ITemporal, which serves as an interface for the temporal
		 *   discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		ITemporal(const CConfig *config_container);
		

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~ITemporal(void);

		/*!
		 * @brief Pure virtual function that computes the upcoming solution in time. Must be overridden
		 *   by a derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 * @param[in] iteration_container input iteration container.
		 * @param[in] openmp_container OpenMP shared parallelization container.
		 * @param[in,out] multizone_solver_container input multizone solver container.
		 * @param[in] physicaltime current physical simulation time.
		 * @param[in] timestep physical time step.
		 */
		virtual void UpdateTime(CConfig            *config_container,
				                    CMultizoneGeometry *multizone_geometry_container,
														CIteration         *iteration_container,
														COpenMP            *openmp_container,
														CMultizoneSolver   *multizone_solver_container,
														as3double           physicaltime, 
														as3double           timestep) = 0;

		/*!
		 * @brief Pure virtual function that updates one element for a temporal stage, without mesh
		 *   traversal.
		 *
		 * @param[in,out] physical_element solution, residual, and temporal history for one element.
		 * @param[in] timestep physical time step.
		 * @param[in] iStageRK index of the RK stage.
		 *
		 * @note Element-local operation with no mesh traversal or synchronization; the caller owns this
		 *   element for the stage.
		 */
		virtual void UpdateElementSolution(CPhysicalElement *physical_element,
		                                     as3double         timestep,
		                                     std::size_t    iStageRK) const = 0;

	protected:

	private:
		// Disable default constructor.
		ITemporal(void) = delete;
		// Disable default copy constructor.
		ITemporal(const ITemporal&) = delete;
		// Disable default copy operator.
		ITemporal& operator=(ITemporal&) = delete;
};
