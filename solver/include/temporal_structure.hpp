#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "iteration_structure.hpp"
#include "solver_structure.hpp"

// Forward declaration to avoid compiler issues.
class ISolver;
class CIteration;

/*!
 * @brief An interface class used for the temporal discretization. 
 */
class ITemporal
{
	public:

		/*!
		 * @brief Constructor of ITemporal, which serves as an interface for the temporal discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		ITemporal(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~ITemporal(void);

		/*!
		 * @brief Pure virtual function that computes the upcoming solution in time. Must be overridden by a derived class.
		 * 
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] iteration_container input iteration container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] timephysical current physical simulation time.
		 * @param[in] timestep physical time step.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		virtual void UpdateTime(CConfig                               *config_container,
				                    CGeometry                             *geometry_container,
														CIteration                            *iteration_container,
														as3vector1d<std::unique_ptr<ISolver>> &solver_container,
														as3double                              physicaltime, 
														as3double                              timestep,
														as3vector1d<as3double>                &monitordata) = 0;

	protected:

	private:
		// Disable default constructor.
		ITemporal(void) = delete;
		// Disable default copy constructor.
		ITemporal(const ITemporal&) = delete;
		// Disable default copy operator.
		ITemporal& operator=(ITemporal&) = delete;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a temporal discretization based on the SSP-RK3. 
 */
class CSSPRK3Temporal final : public ITemporal
{
	public:
	
		/*!
		 * @brief Constructor of CSSPRK3Temporal, which initializes a SSP-RK3 temporal discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CSSPRK3Temporal(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CSSPRK3Temporal(void) final;


		/*!
		 * @brief Function that computes the upcoming solution in time, based on a SSP-RK3. 
		 * 
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] iteration_container input iteration container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] timephysical current physical simulation time.
		 * @param[in] timestep physical time step.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void UpdateTime(CConfig                               *config_container,
		                CGeometry                             *geometry_container,
										CIteration                            *iteration_container,
										as3vector1d<std::unique_ptr<ISolver>> &solver_container,
										as3double                              physicaltime, 
										as3double                              timestep,
										as3vector1d<as3double>                &monitordata) final;

	protected:

	private:

		as3vector1d<as3double> rk3a; ///< SSP-RK3: a-coefficients.
		as3vector1d<as3double> rk3b; ///< SSP-RK3: b-coefficients.
		as3vector1d<as3double> rk3c; ///< SSP-RK3: c-coefficients.


		/*!
		 * @brief Function that evaluates a single stage SSP-RK3. 
		 * 
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] geometry_container input geometry container.
		 * @param[in] iteration_container input iteration container.
		 * @param[in] solver_container input multizone solver container.
		 * @param[in] localtime local time in the RK evaluations.
		 * @param[in] timestep physical time step.
		 * @param[in] alpha relevant SSP-RK3 a-coefficient.
		 * @param[in] beta relevant SSP-RK3 b-coefficient.
		 * @param[out] monitordata vector of parameters to monitor.
		 */
		void EvaluateSSPRK3(CConfig                               *config_container,
		                    CGeometry                             *geometry_container,
										    CIteration                            *iteration_container,
										    as3vector1d<std::unique_ptr<ISolver>> &solver_container,
										    as3double                              localtime, 
										    as3double                              timestep,
												as3double                              alpha,
												as3double                              beta,
												as3vector1d<as3double>                &monitordata);



};


