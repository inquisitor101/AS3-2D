#pragma once

/*!
 * @file
 * @brief Initialization of the conservative Euler solution in each zone.
 */


#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"


// Forward declaration to avoid compiler issues.
class ISolver;

/*!
 * @brief An interface for initializing the conservative solution in one zone.
 */
class IInitialCondition
{
	public:

		/*!
		 * @brief Constructor of IInitialCondition, which serves as an interface for the initial
		 *   condition.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		IInitialCondition(const CConfig *config_container);

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~IInitialCondition(void);


		/*!
		 * @brief A Pure virtual function that initializes the solution over a single zone. Must be
		 *   overridden.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry zone geometry container.
		 * @param[in,out] solver_container solver container.
		 */
		virtual void InitializeSolution(CConfig             *config_container,
				                            CSinglezoneGeometry *singlezone_geometry,
				                            ISolver             *solver_container) = 0;

	protected:

	private:

		// Disable default constructor.
		IInitialCondition(void) = delete;
		// Disable default copy constructor.
		IInitialCondition(const IInitialCondition&) = delete;
		// Disable default copy operator.
		IInitialCondition& operator=(IInitialCondition&) = delete;	
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a Gaussian pressure initial condition.
 */
class CGaussianPressureIC : public IInitialCondition
{
	public:

		/*!
		 * @brief Constructor of CGaussianPressureIC, which initializes a Gaussian pressure IC.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CGaussianPressureIC(const CConfig *config_container);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CGaussianPressureIC(void) override;

		/*!
		 * @brief Function that initializes a Gaussian pressure pulse solution over a single zone.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry zone geometry container.
		 * @param[in,out] solver_container solver container.
		 */
		void InitializeSolution(CConfig             *config_container,
		                        CSinglezoneGeometry *singlezone_geometry,
		                        ISolver             *solver_container) override;

	protected:

	private:
		as3double mXCenter;  ///< Initial disturbance center in the x-direction.
		as3double mYCenter;  ///< Initial disturbance center in the y-direction.

		as3double mRatio;  ///< Initial disturbance strength.
		as3double mWidth;  ///< Initial disturbance width or vortex radius.

		as3double mMachInf;   ///< Freestream Mach number.
		as3double mThetaInf;  ///< Freestream flow angle in degrees.

		as3double mDensityInf;      ///< Freestream density.
		as3double mXVelocityInf;    ///< Freestream velocity in the x-direction.
		as3double mYVelocityInf;    ///< Freestream velocity in the y-direction.
		as3double mPressureInf;     ///< Freestream static pressure.
		as3double mTemperatureInf;  ///< Freestream static temperature.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for an isentropic vortex initial condition.
 */
class CIsentropicVortexIC : public IInitialCondition
{
	public:

		/*!
		 * @brief Constructor of CIsentropicVortexIC, which initializes an isentropic vortex IC.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CIsentropicVortexIC(const CConfig *config_container);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CIsentropicVortexIC(void) override;

		/*!
		 * @brief Function that initializes an isentropic vortex solution over a single zone.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry zone geometry container.
		 * @param[in,out] solver_container solver container.
		 */
		void InitializeSolution(CConfig             *config_container,
		                        CSinglezoneGeometry *singlezone_geometry,
		                        ISolver             *solver_container) override;

	protected:

	private:
		as3double mXCenter;  ///< Initial disturbance center in the x-direction.
		as3double mYCenter;  ///< Initial disturbance center in the y-direction.

		as3double mRatio;  ///< Initial disturbance strength.
		as3double mWidth;  ///< Initial disturbance width or vortex radius.

		as3double mMachInf;   ///< Freestream Mach number.
		as3double mThetaInf;  ///< Freestream flow angle in degrees.

		as3double mDensityInf;      ///< Freestream density.
		as3double mXVelocityInf;    ///< Freestream velocity in the x-direction.
		as3double mYVelocityInf;    ///< Freestream velocity in the y-direction.
		as3double mPressureInf;     ///< Freestream static pressure.
		as3double mTemperatureInf;  ///< Freestream static temperature.
};


