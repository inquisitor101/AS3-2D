#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "input_structure.hpp"


/*!
 * @brief A namespace used for storing specific file import utility functions.
 */
namespace NImportFile
{

	/*!
	 * @brief Function that imports an AS3 grid file.
	 */
	void ImportAS3Grid(CConfig   *config_container, 
			               CGeometry *geometry_container);

	/*!
	 * @brief Function that imports an AS3 grid file in binary format.
	 */
	void ImportAS3GridBinary(CConfig   *config_container, 
			                     CGeometry *geometry_container);

	/*!
	 * @brief Function that checks whether byte-swapping is required.
	 */
	template<typename T>
	bool CheckByteSwapping(const T expect, T test)
	{
		bool swap = ( expect != test ) ? true : false;

		// Check if byte-swapping works or not.
		if( swap )
		{
			T tmp = test;
			NInputUtility::SwapBytes(&tmp, sizeof(T), 1);
			if( expect != tmp )
				ERROR("Issue in file, could not read it.");
		}

		return swap;
	}



}


