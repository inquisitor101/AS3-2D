#pragma once 

#include "config_structure.hpp"



/*!
 * @brief A class for the OpenMP shared memory parallelization. 
 */
class COpenMP
{
	public:

		/*!
		 * @brief Constructor of COpenMP, which initializes the OpenMP class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 */
		COpenMP(CConfig *config_container);
	
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~COpenMP(void);

    size_t GetnThread(void) const { return mNThread; }

	protected:

	private:
    size_t mNThread = 1;

		// Disable default constructor.
		COpenMP(void) = delete;
		// Disable default copy constructor.
		COpenMP(const COpenMP&) = delete;
		// Disable default copy operator.
		COpenMP& operator=(COpenMP&) = delete;
};

