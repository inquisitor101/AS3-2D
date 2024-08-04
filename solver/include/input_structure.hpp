#pragma once

#include "option_structure.hpp"
#include <fstream>
#include <cstring>


/*!
 * @brief A namespace used for storing generic input utility functions.
 */
namespace NInputUtility
{

	/*!
	 * @brief Function that extracts an extension from a filename.
	 *
	 * @param[in] filename input file name.
	 * 
	 * @return extension output extension if detected.
	 */
	std::string GetFileExtension(std::string filename);

  /*!
   * @brief Function that extracts a number after a symbol from a vector of strings.
   *
   * @param[in] str input vector of strings.
   * @param[out] val output value of the integers.
   * @param[in] sym character used as a delimiter.
   */
	template<class TValueType>
	TValueType ExtractLastDigitFromString(std::string  str,
																	      const char  *sym = "_");

  /*!
   * @brief Function that extracts a number after a symbol from a single string.
   *
   * @param[in] str input vector of strings.
   * @param[out] val output value of the integers.
   * @param[in] sym character used as a delimiter.
   */
	template<class TValueType>
	void ExtractLastDigitFromString(as3vector1d<std::string>  str,
			                            as3vector1d<TValueType>  &val,
																	const char               *sym = "_");

	/*!
	 * @brief Function used in reading "scalar" data, takes default parameter.
	 *
	 * @param[in] inputFile reference to input file.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] defaultValue reference to default value specified.
	 * @param[in] ResetLine option to reset to first line of file.
	 */
	template<class TValueType>
	void AddScalarOption(std::ifstream  &inputFile,
	                     const char     *keyword,
	                     TValueType     &value,
	                     TValueType      defaultValue,
	                     bool            ResetLine = false);
	
	/*!
	 * @brief Function used in reading "scalar" data.
	 *
	 * @param[in] inputFile reference to input file.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] ResetLine option to reset to first line of file.
	 */
	template<class TValueType>
	void AddScalarOption(std::ifstream  &inputFile,
	                     const char     *keyword,
	                     TValueType     &value,
	                     bool            ResetLine = false);

	/*!
	 * @brief Function that works as a template to read vector data, otherwise assigns default value.
	 *
	 * @param[in] inputFile reference to input file.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] defaultValue reference to default value specified.
	 * @param[in] ResetLine option to reset to first line of file.
	 */
	template<class TValueType>
	void AddVectorOption(std::ifstream           &inputFile,
	                     const char              *keyword,
	                     std::vector<TValueType> &value,
	                     std::vector<TValueType>  defaultValue,
	                     bool                     ResetLine = false);
	
	/*!
	 * @brief Function that works as a template to read vector data, otherwise exits.
	 *
	 * @param[in] inputFile reference to input file.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] ResetLine option to reset to first line of file.
	 */
	template<class TValueType>
	void AddVectorOption(std::ifstream           &inputFile,
	                     const char              *keyword,
	                     std::vector<TValueType> &value,
	                     bool                     ResetLine = false);

	/*!
	 * @brief Function that swaps bytes.
	 *
	 * @param[in,out] buffer pointer to the data buffer being swapped.
	 * @param[in] nBytes input number of bytes.
	 * @param[in] nItems input number of items in the buffer.
	 */
	void SwapBytes(void   *buffer,
			           size_t  nBytes,
								 size_t  nItems); 

	/*!
	 * @brief Function that removes blanks from a string.
	 *
	 * @param[in] lineBuf reference to input string line.
	 */
	void RemoveBlanks(std::string &lineBuf);
	
	/*!
	 * @brief Function that replaces tabs and return characters in a string with a blank.
	 *
	 * @param[in] lineBuf reference to input string line.
	 */
	void ReplaceTabsAndReturns(std::string &lineBuf);
	
	/*!
	 * @brief Function that returns a lowered-case version of input string.
	 *
	 * @param[in] lineBuf reference to input string line.
	 */
	void CreateLowerCase(std::string &lineBuf);
	
	/*!
	 * @brief Function that searches for and finds a string vector from input file.
	 *
	 * @param[in] inputFile reference to input file.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] valStringVector reference to extracted string vector. 
	 *
	 * @return indicator on whether the extraction worked or not
	 */
	bool FindStringVectorFromKeyword(std::ifstream            &inputFile,
	                          			 const char               *keyword,
	                          			 std::vector<std::string> &valStringVector);
	
	/*!
	 * @brief Function that searches for and finds a string from input file.
	 *
	 * @param[in] inputFile reference to input file.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] valString reference to extracted string. 
	 *
	 * @return indicator on whether the extraction worked or not
	 */
	bool FindStringFromKeyword(std::ifstream &inputFile,
	                           const char     *keyword,
	                           std::string    &valString);

}

// Definitions of the templated functions.
#include "input_structure.inl"


