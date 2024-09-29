#pragma once

#include <iostream>

#ifdef ENABLE_NAN_CHECK
#include <fenv.h>
#endif

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
	/*!
	 * @brief Function that checks floating point errors.
	 */
	void CheckFloatingError(void);
#endif
}
