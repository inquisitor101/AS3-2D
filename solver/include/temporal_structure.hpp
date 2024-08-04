#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"


/*!
 * @brief An interface class used for the temporal discretization. 
 */
class ITemporal
{
	public:
		// Disable default constructor.
		ITemporal(void) = delete;
		
		/*!
		 * @brief Constructor of ITemporal, which serves as an interface for the temporal discretization.
		 */
		ITemporal(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~ITemporal(void);

	protected:

	private:

};

//-----------------------------------------------------------------------------------

/*!
 * @brief A class for a temporal discretization based on the SSP-RK3. 
 */
class CSSPRK3Temporal final : public ITemporal
{
	public:
		// Disable default constructor.
		CSSPRK3Temporal(void) = delete;
		
		/*!
		 * @brief Constructor of CSSPRK3Temporal, which initializes a SSP-RK3 temporal discretization.
		 */
		CSSPRK3Temporal(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CSSPRK3Temporal(void) final;

	protected:

	private:

};


