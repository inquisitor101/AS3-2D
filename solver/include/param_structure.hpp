#pragma once

#include "option_structure.hpp"


/*!
 * @brief A struct used for storing the parameters in boundary markers. 
 */
struct IBoundaryParamMarker
{
	std::string mName; ///< Name of the current marker.
	
	/*!
	 * @brief Virtual destructor, does nothing, but needed for polymorphism.
	 */
	virtual ~IBoundaryParamMarker() {}

	/*!
	 * @brief Pure virtual function that returns the type of BC on this marker. Must be overridden.
	 */
	virtual ETypeBC GetTypeBC(void) const = 0;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A struct used for storing the parameters in periodic markers. 
 */
struct CPeriodicParamMarker : public IBoundaryParamMarker
{
	/*!
	 * @brief Constructor that defines the parameters of this class.
	 * 
	 * @param[in] buffer vector of strings containing the parameters.
	 */
	explicit CPeriodicParamMarker(as3vector1d<std::string> buffer);

	/*!
	 * @brief Function that returns the type of BC on this marker, which is periodic.
	 */
	ETypeBC GetTypeBC(void) const override {return ETypeBC::PERIODIC;}

	std::string              mNameMatching; ///< Name of the matching marker.
  std::array<as3double, 2> mVectorTrans;  ///< Translation vector, from I to J.
};


