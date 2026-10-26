#pragma once 

/*!
 * @file
 * @brief Numerical surface flux interfaces and the Roe Euler flux.
 */


#include "option_structure.hpp"
#include "config_structure.hpp"


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
