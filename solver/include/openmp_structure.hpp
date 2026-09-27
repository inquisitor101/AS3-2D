#pragma once 

#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"

// Forward declaration to avoid compiler issues.
class ISolver;


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
		 * @param[in] geometry_container input geometry container.
		 * @param[in] solver_container input multizone solver container.
		 */
		COpenMP(CConfig                               *config_container,
				    CGeometry                             *geometry_container,
						as3vector1d<std::unique_ptr<ISolver>> &solver_container);
	
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~COpenMP(void);

	protected:

	private:
			
		// Disable default constructor.
		COpenMP(void) = delete;
		// Disable default copy constructor.
		COpenMP(const COpenMP&) = delete;
		// Disable default copy operator.
		COpenMP& operator=(COpenMP&) = delete;
};

