#pragma once

#include <iostream>

/*!
 * @brief Macro for an error log message, before terminating the program.
 */
#define ERROR( msg ) NError::Terminate( __FUNCTION__, __FILE__, __LINE__, msg )



/*!
 * @brief A namespace used for storing error functionalities.
 */
namespace NError
{

	/*!
	 * @brief Function that terminates the program and outputs an error message.
	 *
	 * @param[in] functionName name of the function encountered the error.
	 * @param[in] fileName name of the file which reported the error.
	 * @param[in] lineNumber line number where the error is found.
	 * @param[in] errorMessage log message printed.
	 */
	void Terminate(const char        *functionName,
	               const char        *fileName,
	               const int          lineNumber,
	               const std::string &errorMessage);


#ifdef ENABLE_NAN_CHECK
	void CheckFloatingError
	(
	  void
	)
	 /*
	  * Function that determines if a floating-point error is encountered.
	  */
	{
	  std::string message;
	  // Check what error has been raised.
	  if( std::fetestexcept(FE_DIVBYZERO) ) message += "Error... FE_DIVBYZERO,";
	  // if( std::fetestexcept(FE_INEXACT)   ) message += "Error... FE_INEXACT,";
	  if( std::fetestexcept(FE_INVALID)   ) message += "Error... FE_INVALID,";
	  if( std::fetestexcept(FE_OVERFLOW)  ) message += "Error... FE_OVERFLOW,";
	  //if( std::fetestexcept(FE_UNDERFLOW) ) message += "Error... FE_UNDERFLOW,";
	
	  if( !message.empty() ) ERROR(message);
	}
#endif
}
