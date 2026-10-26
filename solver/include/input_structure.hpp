#pragma once

/*!
 * @file
 * @brief Configuration text parsing, string conversion and byte swapping.
 */


#include <cstddef>
#include <string>
#include <vector>

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
	 * @brief Function that parses the unsigned numeric suffix after the first delimiter in a single input
	 *   string.
	 *
	 * @tparam TValueType type used for parsed or padded values.
	 *
	 * @param[in] str input string.
	 * @param[in] sym one-character delimiter preceding the numeric suffix.
	 *
	 * @return unsigned numeric suffix converted to TValueType; terminates if no delimiter is present.
	 */
	template<class TValueType>
	TValueType ExtractLastDigitFromString(std::string  str,
																	      const char  *sym = "_");

	/*!
	 * @brief Function that parses the unsigned numeric suffix after the first delimiter in each input
	 *   string.
	 *
	 * @tparam TValueType type used for parsed or padded values.
	 *
	 * @param[in] str input strings.
	 * @param[out] val output value of the integers.
	 * @param[in] sym one-character delimiter preceding the numeric suffix.
	 */
	template<class TValueType>
	void ExtractLastDigitFromString(as3vector1d<std::string>  str,
			                            as3vector1d<TValueType>  &val,
																	const char               *sym = "_");

	/*!
	 * @brief Function used in reading "scalar" data, takes default parameter.
	 *
	 * @tparam TValueType type used for parsed or padded values.
	 *
	 * @param[in,out] inputFile input stream whose read position may change.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] defaultValue reference to default value specified.
	 * @param[in] ResetLine whether to restore the original stream position after reading.
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
	 * @tparam TValueType type used for parsed or padded values.
	 *
	 * @param[in,out] inputFile input stream whose read position may change.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] ResetLine whether to restore the original stream position after reading.
	 */
	template<class TValueType>
	void AddScalarOption(std::ifstream  &inputFile,
	                     const char     *keyword,
	                     TValueType     &value,
	                     bool            ResetLine = false);

	/*!
	 * @brief Function that works as a template to read vector data, otherwise assigns default value.
	 *
	 * @tparam TValueType type used for parsed or padded values.
	 *
	 * @param[in,out] inputFile input stream whose read position may change.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] defaultValue reference to default value specified.
	 * @param[in] ResetLine whether to restore the original stream position after reading.
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
	 * @tparam TValueType type used for parsed or padded values.
	 *
	 * @param[in,out] inputFile input stream whose read position may change.
	 * @param[in] keyword pointer to keyword looked for.
	 * @param[out] value reference to value read.
	 * @param[in] ResetLine whether to restore the original stream position after reading.
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
			           std::size_t  nBytes,
								 std::size_t  nItems); 

	/*!
	 * @brief Function that removes blanks from a string.
	 *
	 * @param[in,out] lineBuf reference to input string line.
	 */
	void RemoveBlanks(std::string &lineBuf);
	

	/*!
	 * @brief Function that replaces tabs and return characters in a string with a blank.
	 *
	 * @param[in,out] lineBuf reference to input string line.
	 */
	void ReplaceTabsAndReturns(std::string &lineBuf);
	

	/*!
	 * @brief Function that converts the input string to lower case in place.
	 *
	 * @param[in,out] lineBuf reference to input string line.
	 */
	void CreateLowerCase(std::string &lineBuf);
	

	/*!
	 * @brief Function that searches for and finds a string vector from input file.
	 *
	 * @param[in,out] inputFile input stream whose read position may change.
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
	 * @param[in,out] inputFile input stream whose read position may change.
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


