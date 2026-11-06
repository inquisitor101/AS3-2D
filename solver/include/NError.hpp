#pragma once

/*!
 * @file
 * @brief Error reporting, warnings and optional floating-point exception checks.
 */


#include <string>

#include <iostream>
#include <cstdlib>

#ifdef ENABLE_NAN_CHECK
#include <cfenv>
#endif

/*!
 * @brief Macro for an error log message, before terminating the program.
 *
 * @param[in] msg diagnostic message printed with its source location.
 */
#define ERROR( msg ) NError::Terminate( __FUNCTION__, __FILE__, __LINE__, msg )

/*!
 * @brief Macro that reports a warning and continues execution.
 *
 * @param[in] msg diagnostic message printed with its source location.
 */
#define WARNING( msg ) NError::Warning( __FUNCTION__, __FILE__, __LINE__, msg )




/*!
 * @brief A namespace used for storing error functionalities.
 */
namespace NError
{

	/*!
	 * @brief Function that terminates the program and outputs an error message.
	 *
	 * @param[in] func name of the function encountered the error.
	 * @param[in] file name of the file which reported the error.
	 * @param[in] line line number where the error is found.
	 * @param[in] error log message printed.
	 */
	[[noreturn]] 
	void Terminate(const char        *func,
	               const char        *file,
	               const int          line,
	               const std::string &error);

	/*!
	 * @brief Function that issues a warning message in the program.
	 *
	 * @param[in] func name of the function encountered the warning.
	 * @param[in] file name of the file which reported the warning.
	 * @param[in] line line number where the warning is found.
	 * @param[in] warning log message printed.
	 */
	void Warning(const char        *func,
	             const char        *file,
	             const int          line,
	             const std::string &warning);


#ifdef ENABLE_NAN_CHECK

	/*!
	 * @brief Preserve solver exception flags across external numerical library calls.
	 */
	class CLibraryFloatingPointScope
	{
	  public:
	  	CLibraryFloatingPointScope();
	  	~CLibraryFloatingPointScope();
	  	CLibraryFloatingPointScope(const CLibraryFloatingPointScope&) = delete;
	  	CLibraryFloatingPointScope& operator=(const CLibraryFloatingPointScope&) = delete;

	  private:
	  	std::fenv_t mEnvironment;
	};

	/*!
	 * @brief Function that checks floating point errors.
	 */
	void CheckFloatingError(void);
#endif
}
