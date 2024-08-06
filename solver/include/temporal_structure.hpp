#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"


/*!
 * @brief An interface class used for the temporal discretization. 
 */
class ITemporal
{
	public:

		/*!
		 * @brief Constructor of ITemporal, which serves as an interface for the temporal discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		ITemporal(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~ITemporal(void);

	protected:

	private:
		// Disable default constructor.
		ITemporal(void) = delete;
		// Disable default copy constructor.
		ITemporal(const ITemporal&) = delete;
		// Disable default copy operator.
		ITemporal& operator=(ITemporal&) = delete;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a temporal discretization based on the SSP-RK3. 
 */
class CSSPRK3Temporal final : public ITemporal
{
	public:
	
		/*!
		 * @brief Constructor of CSSPRK3Temporal, which initializes a SSP-RK3 temporal discretization.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CSSPRK3Temporal(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CSSPRK3Temporal(void) final;

	protected:

	private:
		as3vector1d<as3double> rk3a; ///< SSP-RK3: a-coefficients.
		as3vector1d<as3double> rk3b; ///< SSP-RK3: b-coefficients.
		as3vector1d<as3double> rk3c; ///< SSP-RK3: c-coefficients.


};


