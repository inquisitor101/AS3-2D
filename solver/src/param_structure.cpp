#include "param_structure.hpp" 


//-----------------------------------------------------------------------------------
// CInterfaceParamMarker member functions.
//-----------------------------------------------------------------------------------


CInterfaceParamMarker::CInterfaceParamMarker
(
 as3vector1d<std::string> buffer
)
 /*
	* Constructor that initializes the parameters needed in an interface BC.
	*/
{
	// Deduce the translation vector from the buffer string and store it into tmp.
	as3double tmp[2];
	std::istringstream istr( buffer[2] + " " + buffer[3] );
	istr >> tmp[0] >> tmp[1];
	
	// Check if any bad character is input.
	if( istr.fail() ) ERROR("Wrong character input inside translation vector.");

	// This marker is from i's perspective: i to j.
	mName           = buffer[0];
	mNameMatching   = buffer[1];
	mVectorTrans[0] = tmp[0];
	mVectorTrans[1] = tmp[1];
	
	// Ensure the boundary markers are unique.
	if( mName == mNameMatching ) ERROR("Interface markers must be unique.");
}


