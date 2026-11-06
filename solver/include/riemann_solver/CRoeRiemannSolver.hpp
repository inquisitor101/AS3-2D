#pragma once

/*!
 * @file
 * @brief Roe numerical flux for the Euler equations.
 */


#include "riemann_solver/IRiemannSolver.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"


/*!
 * @brief A class for a Riemann solver based on Roe's method.
 */
class CRoeRiemannSolver final: public IRiemannSolver
{
	public:
	

		/*!
		 * @brief Constructor of CRoeRiemannSolver, which initializes Roe's Riemann solver.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CRoeRiemannSolver(const CConfig *config_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CRoeRiemannSolver(void) override;

		/*!
		 * @brief Getter function which returns the numerical flux method.
		 *
		 * @return implemented numerical flux method.
		 */
		ETypeRiemannSolver GetTypeRiemannSolver(void) const override {return ETypeRiemannSolver::ROE;}

		/*!
		 * @brief Function that computes the weighted Euler numerical surface flux using Roe's method.
		 *
		 * @param[in] wts integration weights on the standard element.
		 * @param[in] met surface length factor and supplied unit normal, stored as [3][point].
		 * @param[in] solL left conservative state [4][point]: density, x-momentum, y-momentum, total
		 *   energy density.
		 * @param[in] solR right conservative state [4][point] in the same variable order as solL.
		 * @param[out] flux numerical flux [4][point], including quadrature weight, surface length and
		 *   residual sign.
		 */
		void ComputeFlux(const CMatrixAS3<as3double>     &wts,
		                 const CMatrixAS3<as3double>     &met,
						         const CWorkMatrixAS3<as3double> &solL,
		                 const CWorkMatrixAS3<as3double> &solR,
						         CWorkMatrixAS3<as3double>       &flux) override;

	protected:

	private:
		const as3double mDelta = 0.001;  ///< Entropy fix constant.
};
