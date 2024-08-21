#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "solver_structure.hpp"


/*!
 * @brief A namespace used for storing logging functionalities.
 */
namespace NLogger
{

	/*!
	 * @brief Function that prints the entries of a matrix.
	 *
	 * @param[in] m reference to the matrix.
	 * @param[in] g log message printed.
	 */
	template<typename T>
	void PrintMatrixAS3(CMatrixAS3<T> &m, std::string g)
	{
		std::cout << "\n" << g << "\n";
		for(size_t i=0; i<m.row(); i++)
		{
			for(size_t j=0; j<m.col(); j++) std::cout << m(i,j) << " ";
			std::cout << std::endl;
		}
	}

	/*!
	 * @brief Function that prints the type of solver in each zone.
	 *
	 * @param[in] config_container configuration/dictionary container.
	 */
	void PrintInitSolver(CConfig *config_container);

	/*!
	 * @brief Function that displays the boundary condition information over all zones.
	 *
	 * @param[in] config_container configuration/dictionary container.
	 * @param[in] solver_container input vector of solver containers.
	 */
	void DisplayBoundaryConditions(CConfig                               *config_container,
			                           as3vector1d<std::unique_ptr<ISolver>> &solver_container);
}
