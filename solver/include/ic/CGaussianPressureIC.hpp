#pragma once

/*!
 * @file
 * @brief Gaussian pressure pulse initial condition.
 */


#include "common/AS3ScalarTypes.hpp"
#include "ic/IInitialCondition.hpp"


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
