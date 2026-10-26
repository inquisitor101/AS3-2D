#pragma once 

/*!
 * @file
 * @brief Monitoring values collected during simulation advancement.
 */


#include <cstddef>

#include "option_structure.hpp"
#include "config_structure.hpp"


/*!
 * @brief Monitoring values accumulated during time-step selection and advancement.
 */
struct CMonitorData
{

	/*!
	 * @brief Constructor that defines this class.
	 *
	 * @param[in] config_container configuration/dictionary container.
	 */
	CMonitorData(const CConfig *config_container);

	/*!
	 * @brief Destructor that releases the resources owned by this object.
	 */
	~CMonitorData(void) = default;

	as3double mMachMax      = C_ZERO;  ///< Max Mach number.
	as3double mDensityL2    = C_ZERO;  ///< Reserved density L2 diagnostic; not currently updated.
	as3double mMinTimeStep  = C_ZERO;  ///< Minimum time step.
	as3double mMaxTimeStep  = C_ZERO;  ///< Maximum time step.
	std::size_t    mNSyncSubStep = 0;  ///< Number of synchronization sub-steps.
};


