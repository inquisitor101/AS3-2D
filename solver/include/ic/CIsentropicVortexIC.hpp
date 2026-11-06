#pragma once

/*!
 * @file
 * @brief Isentropic vortex initial condition.
 */


#include "common/AS3ScalarTypes.hpp"
#include "ic/IInitialCondition.hpp"


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
		                        ISinglezoneSolver  *solver_container) override;

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
