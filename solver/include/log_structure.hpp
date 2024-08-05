#pragma once


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

}