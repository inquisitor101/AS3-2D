#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"


/*!
 * @brief A class for a single grid sweep across all zones. 
 */
class CIteration
{
	public:

		/*!
		 * @brief Constructor of CIteration, which is responsible for an iteration via a single grid sweep.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		CIteration(CConfig *config_container);
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CIteration(void);

		/*!
		 * @brief Function that performs a grid sweep over all the zones.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		void GridSweep(CConfig *config_container);

	protected:

	private:

		/*!
		 * @brief Function that preprocesses the solution, before sweeping the grid.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		void PreProcessIteration(CConfig *config_container);

		/*!
		 * @brief Function that postprocesses the solution, after sweeping the grid.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		void PostProcessIteration(CConfig *config_container);


		// Disable default constructor.
		CIteration(void) = delete;
		// Disable default copy constructor.
		CIteration(const CIteration&) = delete;
		// Disable default copy operator.
		CIteration& operator=(CIteration&) = delete;
};


