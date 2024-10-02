#include "error_structure.hpp"


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
  * Function which prints an error message and exit from the program.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the error message
  // looks nicer.
  std::string  message = error;
  std::string  buffer  = message;
  std::string::size_type off = 1;

  for(;;)
  {
    std::string::size_type loc = buffer.find("\n");
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
  exit(1);
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
  * Function which prints a warning message, before resuming the program.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the warning message
  // looks nicer.
  std::string  message = warning;
  std::string  buffer  = message;
  std::string::size_type off = 1;

  for(;;)
  {
    std::string::size_type loc = buffer.find("\n");
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
void NError::CheckFloatingError
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

