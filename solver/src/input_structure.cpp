#include "input_structure.hpp"

#include <cctype>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <vector>



//-----------------------------------------------------------------------------------
// NInputUtility namespace functions.
//-----------------------------------------------------------------------------------

std::string NInputUtility::GetFileExtension
(
 std::string filename
)
 /*
	* Function that extracts the filename extension.
	*/
{
	// Check whether or not an extension is present.
	if( filename.find_last_of(".") != std::string::npos )
	{
		// Extension is detected.
		return filename.substr( filename.find_last_of(".")+1 );
	}
	else
	{
		// No extension is found.
		return "";
	}
}

//-----------------------------------------------------------------------------------

void NInputUtility::SwapBytes
(
 void   *buffer,
 std::size_t  nBytes,
 std::size_t  nItems
)
 /*
	* Function that reverses the byte order of each item in a buffer.
	*/
{
	// Store half the number of bytes in kk and cast the buffer
	// to a character buffer.
	char *buf       = (char *) buffer;
	const std::size_t kk = nBytes/2;
	
	// Loop over the number of items in the buffer.
	for(std::size_t j=0; j<nItems; ++j)
	{
	  // Initialize ii and jj, which are used to store the
	  // indices of the bytes to be swapped.
	  std::size_t ii = j*nBytes;
	  std::size_t jj = ii + nBytes - 1;
	
	  // Swap the bytes.
	  for(std::size_t i=0; i<kk; ++i, ++ii, --jj)
	  {
	    const char tmp = buf[jj];
	    buf[jj] = buf[ii];
	    buf[ii] = tmp;
	  }
	}
}

//-----------------------------------------------------------------------------------

void NInputUtility::RemoveBlanks
(
 std::string &lineBuf
)
 /*
  * Function that removes spaces from a string in place.
  */
{
  std::size_t first = 0;
  std::size_t last = lineBuf.size();
  while(first < last && std::isspace(static_cast<unsigned char>(lineBuf[first]))) ++first;
  while(last > first && std::isspace(static_cast<unsigned char>(lineBuf[last-1]))) --last;
  lineBuf = lineBuf.substr(first, last-first);
}

//-----------------------------------------------------------------------------------

void NInputUtility::ReplaceTabsAndReturns
(
 std::string &lineBuf
)
 /*
  * Function that replaces tabs and line endings with spaces in place.
  */
{
  // Replace the tabs.
  for(;;)
  {
    std::size_t pos = lineBuf.find("\t");
    if(pos == std::string::npos) break;
    lineBuf.replace(pos, 1, " ");
  }

  // Replace the returns.
  for(;;)
  {
    std::size_t pos = lineBuf.find("\n");
    if(pos == std::string::npos) break;
    lineBuf.replace(pos, 1, " ");
  }
}


void NInputUtility::CreateLowerCase
(
 std::string &lineBuf
)
 /*
  * Function that converts a string to lower case in place.
  */
{
  // Determine the length of the string and convert its elements to
  // lower case.
  std::size_t strLength = lineBuf.length();
  for(std::size_t i=0; i<strLength; ++i)
    lineBuf[i] = (char) tolower(lineBuf[i]);
}

//-----------------------------------------------------------------------------------

bool NInputUtility::FindStringVectorFromKeyword
(
 std::ifstream            &inputFile,
 const char               *keyword,
 std::vector<std::string> &valStringVector
)
/*
 * Function that finds a configuration keyword and extracts its vector entries.
 */
{
  // Length of target string.
  const std::size_t lenKeyword = std::strlen(keyword);

  // Initialize search as false.
  bool keywordFound = false;

  // Line stream.
  std::string lineBuf;

  // Read line by line, until/if value of string is matched.
  while(std::getline(inputFile, lineBuf)){
    // Remove tabs and new lines.
    ReplaceTabsAndReturns(lineBuf);

    // Ignore comments, i.e. "%"
    std::size_t pos = lineBuf.find("%");
    if(pos != std::string::npos ) lineBuf.erase(pos, lineBuf.length()-pos);

    // Remove any blanks in line.
    RemoveBlanks(lineBuf);

    // Evaluate current line.
    if( lineBuf.length() ){

      // Compare line with target string.
      if(lineBuf.compare(0, lenKeyword, keyword) == 0){
				// Find symbol ( position.
				pos = lineBuf.find("(");

				// Current line found begins after ( symbol.
				std::string lineFound = lineBuf.substr(pos+1);

				// Initialize string stream of current line.
				std::stringstream lineStream(lineFound);
				// Iterate while the current line is not over.
				while( lineStream.good() ) {

					// Initialize data.
					std::string data;
					// Populate data.
					std::getline(lineStream, data, ',');
					// Remove blanks in keyword.
					RemoveBlanks(data);
					// Add found value to vector.
					valStringVector.push_back(data);

					// Check whether data is corrupted or empty.
					if( data == "" ) ERROR("Data is corrupted and/or empty!");
				}

				// If string not empty, assert the format used is appropriate.
				if( !valStringVector.empty() ){

					// Extract final index.
					const std::size_t idx = valStringVector.size()-1;

					// Extract last character of string.
					char endSymbol = valStringVector[idx].back();

					// Exit immediately, if incosistent use is found.
					if( endSymbol != ')' ) ERROR("Wrong format used, proper use is: ( data1, data2, etc )");

					// Remove final character: right-parenthesis symbol.
					valStringVector[idx].pop_back();
				}

        // Keyword found.
        keywordFound = true;
      }
    }
    // If target string is found, break out of loop.
    if(keywordFound) break;
  }

  // Return happily.
  return keywordFound;
}

//-----------------------------------------------------------------------------------

bool NInputUtility::FindStringFromKeyword
(
 std::ifstream &inputFile,
 const char    *keyword,
 std::string   &valString
)
/*
 * Function that finds a configuration keyword and extracts its scalar text.
 */
{
  // Length of target string.
  const std::size_t lenKeyword = std::strlen(keyword);

  // Initialize search as false.
  bool keywordFound = false;

  // Line stream.
  std::string lineBuf;

  // Read line by line, until/if value of string is matched.
  while(std::getline(inputFile, lineBuf)){
    // Remove tabs and new lines.
    ReplaceTabsAndReturns(lineBuf);

    // Ignore comments, i.e. "%"
    std::size_t pos = lineBuf.find("%");
    if(pos != std::string::npos ) lineBuf.erase(pos, lineBuf.length()-pos);

    // Remove any blanks in line.
    RemoveBlanks(lineBuf);

    // Evaluate current line.
    if( lineBuf.length() ){
      // Create a lower-case version of current line.
      std::string lineBufLower = lineBuf;
      //CreateLowerCase(lineBufLower);

      // Compare line with target string.
      if(lineBufLower.compare(0, lenKeyword, keyword) == 0){
        pos = lineBuf.find("=");
        valString = lineBuf.substr(pos+1);

        // Keyword found.
        keywordFound = true;
      }
    }
    // If target string is found, break out of loop.
    if(keywordFound) break;
  }

  // Return happily.
  return keywordFound;
}


