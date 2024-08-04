#include "error_log.hpp"



void Terminate
(
  const char        *functionName,
  const char        *fileName,
  const int          lineNumber,
  const std::string &errorMessage
)
 /*
  * Function which prints an error message and exit from the program.
  */
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the error message
  // looks nicer.
  std::string  message = errorMessage;
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
	std::cout << std::endl;
	std::cout << "#" << std::endl
            << "#============================= !!! Error !!! "
            << "==============================" << std::endl
            << "#" << std::endl;

  // Write the function name, file name and line number.
  std::cout << "#* Run-time error in " << functionName << std::endl;
  std::cout << "#* File: " << fileName
            << ", Line: "  << lineNumber << std::endl
            << "#" << std::endl;

  // Write the error message and the terminating line.
  std::cout << "# " << message << std::endl
            << "#" << std::endl;
  std::cout << "#============================================"
            << "==============================" << std::endl
            << "#" << std::endl << std::flush;

  // And exit.
  exit(1);
}


