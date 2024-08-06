#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"


// Forward declaration to avoid compiler issues.
class ISolver;

/*!
 * @brief An interface class used for storing initializing the solution.
 */
class IInitialCondition
{
	public:

		/*!
		 * @brief Constructor of IInitialCondition, which serves as an interface for the initial condition.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		IInitialCondition(CConfig *config_container);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IInitialCondition(void);


		/*!
		 * @brief A Pure virtual function that initializes the solution over a single zone. Must be overridden.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] zone_geometry zone geometry container. 
		 * @param[in] solver_container solver container.
		 */
		virtual void InitializeSolution(CConfig       *config_container,
				                            CZoneGeometry *zone_geometry,
				                            ISolver       *solver_container) = 0;

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
		CGaussianPressureIC(CConfig *config_container);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CGaussianPressureIC(void) override;

		/*!
		 * @brief Function that initializes a Gaussian pressure pulse solution over a single zone.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] zone_geometry zone geometry container. 
		 * @param[in] solver_container solver container.
		 */
		void InitializeSolution(CConfig       *config_container,
		                        CZoneGeometry *zone_geometry,
		                        ISolver       *solver_container) override;

	protected:

	private:
		as3double mXCenter;
		as3double mYCenter;

		as3double mRatio; 
		as3double mWidth;

		as3double mMachInf;
		as3double mThetaInf;

		as3double mDensityInf;
		as3double mXVelocityInf;
		as3double mYVelocityInf;
		as3double mPressureInf;
		as3double mTemperatureInf;
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
		CIsentropicVortexIC(CConfig *config_container);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CIsentropicVortexIC(void) override;

		/*!
		 * @brief Function that initializes an isentropic vortex solution over a single zone.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] zone_geometry zone geometry container. 
		 * @param[in] solver_container solver container.
		 */
		void InitializeSolution(CConfig       *config_container,
		                        CZoneGeometry *zone_geometry,
		                        ISolver       *solver_container) override;

	protected:

	private:
		as3double mXCenter;
		as3double mYCenter;

		as3double mRatio; 
		as3double mWidth;

		as3double mMachInf;
		as3double mThetaInf;

		as3double mDensityInf;
		as3double mXVelocityInf;
		as3double mYVelocityInf;
		as3double mPressureInf;
		as3double mTemperatureInf;
};


