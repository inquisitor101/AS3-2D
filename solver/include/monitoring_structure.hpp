#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"


/*!
 * @brief A struct used for storing the monitoring data. 
 */
struct CMonitorData
{
	/*!
	 * @brief Constructor that defines this class.
	 * 
	 * @param[in] config_container configuration/dictionary container.
	 */
	CMonitorData(CConfig *config_container);

	// Default destructor.
	~CMonitorData(void) = default;

	as3double mMachMax      = C_ZERO; ///< Max Mach number.
	as3double mDensityL2    = C_ZERO; ///< L2 norm of the density.
	as3double mMinTimeStep  = C_ZERO; ///< Minimum time step.
	as3double mMaxTimeStep  = C_ZERO; ///< Maximum time step.
	size_t    mNSyncSubStep = 0;      ///< Number of synchronization sub-steps.
};


