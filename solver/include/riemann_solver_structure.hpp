#pragma once 

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
		IRiemannSolver(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IRiemannSolver(void);

		/*!
		 * @brief Pure virtual function that computes a unique flux state on a surface. Must be overridden.
		 *
		 * @param[in] wts integration weights on the standard element.
		 * @param[in] met metrics on this face, based on the owner element/side.
		 * @param[in] solL solution on the left face.
		 * @param[in] solR solution on the right face.
		 * @param[out] flux unique flux at this face.
		 */
		virtual void ComputeFlux(const CMatrixAS3<as3double>     &wts,
				                     const CMatrixAS3<as3double>     &met,
										         const CWorkMatrixAS3<as3double> &solL,
				                     const CWorkMatrixAS3<as3double> &solR,
										         CWorkMatrixAS3<as3double>       &flux) = 0;

		/*!
		 * @brief Pure virtual function that returns the type of Riemann solver. Must be overridden.
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
		CRoeRiemannSolver(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CRoeRiemannSolver(void) final;

		/*!
		 * @brief Function that returns the type of Riemann solver.
		 */
		ETypeRiemannSolver GetTypeRiemannSolver(void) const {return ETypeRiemannSolver::ROE;}

		/*!
		 * @brief Function that computes a unique flux state on a surface using Roe's method.
		 *
		 * @param[in] wts integration weights on the standard element.
		 * @param[in] met metrics on this face, based on the owner element/side.
		 * @param[in] solL solution on the left face.
		 * @param[in] solR solution on the right face.
		 * @param[out] flux unique flux at this face.
		 */
		void ComputeFlux(const CMatrixAS3<as3double>     &wts,
		                 const CMatrixAS3<as3double>     &met,
						         const CWorkMatrixAS3<as3double> &solL,
		                 const CWorkMatrixAS3<as3double> &solR,
						         CWorkMatrixAS3<as3double>       &flux) final;

	protected:

	private:
		const as3double mDelta = 0.001;  ///< Entropy fix constant.
};
