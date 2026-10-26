#include "error_structure.hpp"

#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <ostream>
#include <string>


//-----------------------------------------------------------------------------------
// NError namespace functions.
//-----------------------------------------------------------------------------------


void NError::Terminate
(
 const char        *func,
 const char        *file,
 const int          line,
 const std::string &error
)
 /*
  * Function that prints a source-located error message and terminates execution.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the error message
  // looks nicer.
  std::string  message = error;
  std::string  buffer  = message;
  std::size_t off = 1;

  for(;;)
  {
    std::size_t loc = buffer.find("\n");
    if(loc == std::string::npos) break;
    message.insert(loc+off, "# ");
    off += loc+3;
    buffer.erase(0,loc+1);
  }

  // Header of the error message.
	std::cout << "\n#\n"
            << "#============================= !!! Error !!! "
            << "==============================\n"
            << "#" << std::endl;

  // Write the function name, file name and line number.
  std::cout << "# Run-time error in...\n";
  std::cout << "#   (*) Function: " << func << "\n"
		        << "#   (*) Filename: " << file << "\n"
            << "#   (*) Line num: " << line << "\n"
            << "#"                  << std::endl;

  // Write the warning error and the terminating line.
  std::cout << "#   (!) Reason: \n# " << message << "\n#\n";
  std::cout << "#============================================"
            << "==============================\n"
            << "#"  << std::endl;

  // And exit.
  exit(EXIT_FAILURE);
}

//-----------------------------------------------------------------------------------

void NError::Warning
(
 const char        *func,
 const char        *file,
 const int          line,
 const std::string &warning
)
 /*
  * Function that prints a source-located warning and continues execution.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the warning message
  // looks nicer.
  std::string  message = warning;
  std::string  buffer  = message;
  std::size_t off = 1;

  for(;;)
  {
    std::size_t loc = buffer.find("\n");
    if(loc == std::string::npos) break;
    message.insert(loc+off, "# ");
    off += loc+3;
    buffer.erase(0,loc+1);
  }

  // Header of the warning message.
	std::cout << "\n#\n"
            << "#============================ !!! Warning !!! "
            << "=============================\n"
            << "#" << std::endl;

  // Write the function name, file name and line number.
  std::cout << "# Warning encountered in...\n";
  std::cout << "#   (*) Function: " << func << "\n"
		        << "#   (*) Filename: " << file << "\n"
            << "#   (*) Line num: " << line << "\n"
            << "#"              << std::endl;

  // Write the warning message and the terminating line.
  std::cout << "#   (!) Reason: \n# " << message << "\n#\n";
  std::cout << "#============================================"
            << "==============================\n"
            << "#"  << std::endl;
}

//-----------------------------------------------------------------------------------

#ifdef ENABLE_NAN_CHECK

NError::CLibraryFloatingPointScope::CLibraryFloatingPointScope
(
 void
)
 /*
  * Constructor that saves the caller's floating-point environment before a library call.
  */
{
  // Save the current environment and clear exception flags for the library call.
  if( std::feholdexcept(&mEnvironment) != 0 ) ERROR("Could not save the floating-point environment.");
}

//-----------------------------------------------------------------------------------

NError::CLibraryFloatingPointScope::~CLibraryFloatingPointScope
(
 void
)
 /*
  * Destructor that restores the caller's floating-point environment after a library call.
  */
{
  // Restore the saved environment, including any exception flags raised by the caller.
  if( std::fesetenv(&mEnvironment) != 0 ) ERROR("Could not restore the floating-point environment.");
}

//-----------------------------------------------------------------------------------

void NError::CheckFloatingError
(
  void
)
 /*
  * Function that checks and reports floating-point exception flags.
  */
{
  std::string message;
  // Check what error has been raised.
  if( std::fetestexcept(FE_DIVBYZERO) ) message += "Error... FE_DIVBYZERO,";
  // if( fetestexcept(FE_INEXACT)   ) message += "Error... FE_INEXACT,";
  if( std::fetestexcept(FE_INVALID)   ) message += "Error... FE_INVALID,";
  if( std::fetestexcept(FE_OVERFLOW)  ) message += "Error... FE_OVERFLOW,";
  //if( fetestexcept(FE_UNDERFLOW) ) message += "Error... FE_UNDERFLOW,";

  if( !message.empty() ) ERROR(message);

  // Clear the checked flags before the next solver operation.
  std::feclearexcept(FE_ALL_EXCEPT);
}
#endif

