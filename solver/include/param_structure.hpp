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
 * @brief A struct used for storing the parameters in interface/periodic markers. 
 */
struct CInterfaceParamMarker : public IBoundaryParamMarker
{
	/*!
	 * @brief Constructor that defines the parameters of this class.
	 * 
	 * @param[in] buffer vector of strings containing the parameters.
	 */
	explicit CInterfaceParamMarker(as3vector1d<std::string> buffer);

	/*!
	 * @brief Function that returns the type of BC on this marker, which is an interface.
	 */
	ETypeBC GetTypeBC(void) const override {return ETypeBC::INTERFACE;}

	std::string mNameMatching;   ///< Name of the matching marker.
	as3double   mVectorTrans[2]; ///< Translation vector, from I to J.
};


