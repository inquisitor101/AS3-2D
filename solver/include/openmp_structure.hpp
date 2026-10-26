#pragma once 

/*!
 * @file
 * @brief Configuration and reporting of the OpenMP team size.
 */


#include <cstddef>

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
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~COpenMP(void);

    /*!
     * @brief Getter function which returns the configured OpenMP team size.
     *
     * @return configured OpenMP team size.
     */
    std::size_t GetnThread(void) const { return mNThread; }

	protected:

	private:
    std::size_t mNThread = 1; ///< Configured number of OpenMP threads; one in a serial build.

		// Disable default constructor.
		COpenMP(void) = delete;
		// Disable default copy constructor.
		COpenMP(const COpenMP&) = delete;
		// Disable default copy operator.
		COpenMP& operator=(COpenMP&) = delete;
};

