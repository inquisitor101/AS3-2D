#pragma once

/*!
 * @file
 * @brief Numerical surface flux interface.
 */


#include "common/AS3Options.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"

class CConfig;


/*!
 * @brief An interface class used for the Riemann solver specification.
 */
class IRiemannSolver
{
	public:
	

		/*!
		 * @brief Constructor of IRiemannSolver, which serves as an interface for the Riemann solver.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		IRiemannSolver(const CConfig *config_container);
		

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~IRiemannSolver(void);

		/*!
		 * @brief Pure virtual function that computes a weighted numerical surface flux.
		 *
		 * @param[in] wts integration weights on the standard element.
		 * @param[in] met surface length factor and supplied unit normal, stored as [3][point].
		 * @param[in] solL left conservative state [4][point]: density, x-momentum, y-momentum, total
		 *   energy density.
		 * @param[in] solR right conservative state [4][point] in the same variable order as solL.
		 * @param[out] flux numerical flux [4][point], including quadrature weight, surface length and
		 *   residual sign.
		 */
		virtual void ComputeFlux(const CMatrixAS3<as3double>     &wts,
				                     const CMatrixAS3<as3double>     &met,
										         const CWorkMatrixAS3<as3double> &solL,
				                     const CWorkMatrixAS3<as3double> &solR,
										         CWorkMatrixAS3<as3double>       &flux) = 0;

		/*!
		 * @brief Getter function which returns the numerical flux method.
		 *
		 * @return implemented numerical flux method.
		 */
		virtual ETypeRiemannSolver GetTypeRiemannSolver(void) const = 0;

	protected:

	private:
	
		// Disable default constructor.
		IRiemannSolver(void) = delete;
		// Disable default copy constructor.
		IRiemannSolver(const IRiemannSolver&) = delete;
		// Disable default copy operator.
		IRiemannSolver& operator=(IRiemannSolver&) = delete;
};
