/*!
 * @file
 * @brief Inline implementations for configuration text parsing, string conversion and byte swapping.
 */

#include <cstddef>
#include <fstream>
#include <ios>
#include <sstream>
#include <string>
#include <vector>



//-----------------------------------------------------------------------------------
// Implementation of the templated functions in NInputUtility.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that parses the unsigned numeric suffix after the first delimiter in each input string.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in] str string or strings containing a numeric suffix.
 * @param[out] val numeric suffixes extracted from the strings.
 * @param[in] sym delimiter preceding the numeric suffix.
 */
template<class TValueType>
void NInputUtility::ExtractLastDigitFromString
(
 as3vector1d<std::string>  str,
 as3vector1d<TValueType>  &val,
 const char               *sym
)
 /*
	* Function that extracts a number after a symbol from a string.
	*/
{
	// First, resize the data according to input string dimension.
	val.resize( str.size() );

	// Loop over each entry in the string and extract the final number.
	for(std::size_t i=0; i<str.size(); i++)
	{
		val[i] = ExtractLastDigitFromString<TValueType>(str[i], sym);
	}
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that parses the unsigned numeric suffix after the first delimiter in a single input string.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in] str string or strings containing a numeric suffix.
 * @param[in] sym delimiter preceding the numeric suffix.
 *
 * @return unsigned numeric suffix converted to TValueType; terminates if no delimiter is present.
 */
template<class TValueType>
TValueType NInputUtility::ExtractLastDigitFromString
(
 std::string  str,
 const char  *sym
)
 /*
	* Function that extracts a number after a symbol from a string.
	*/
{	
	// Location of the symbol.
	std::size_t pos = str.find(sym);
	
	// If located, extract the number.
	if( pos != std::string::npos )
	{
		return static_cast<TValueType>( std::stoul( str.substr(pos+1) ) );
	}
	else
	{
		// Otherwise, report an error.
		ERROR("The data must be separated by a symbol, e.g. DATA_#");
		
		// To avoid a compiler warning.
		return 0;
	}
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function used in reading "scalar" data, takes default parameter.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in,out] inputFile input stream whose read position may change.
 * @param[in] keyword configuration keyword to find.
 * @param[in] value parsed output value.
 * @param[in] defaultValue value used when the keyword is absent.
 * @param[in] ResetLine whether to restore the original stream position after reading.
 */
template<class TValueType>
void NInputUtility::AddScalarOption
(
 std::ifstream  &inputFile,
 const char     *keyword,
 TValueType     &value,
 TValueType      defaultValue,
 bool            ResetLine
)
 /*
  * Function that assigns a value from input file, otherwise uses default.
  */
{
  // Initialize string.
  std::string valString;
  // Initialize line number.
  std::streampos CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringFromKeyword(inputFile, keyword, valString) )
	{
		// Assign default value.
    value = defaultValue;
    // Go back to previous line location.
    if( ResetLine ) { inputFile.clear(); inputFile.seekg(CurrentLine); }
    return;
  }

  // Assign string value.
  std::istringstream istr(valString);
  istr >> value;

	// Check if any bad character is input.
	if( istr.fail() ) ERROR("Wrong character input inside expected parameter type!");

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}
	
//-----------------------------------------------------------------------------------

/*!
 * @brief Function used in reading "scalar" data.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in,out] inputFile input stream whose read position may change.
 * @param[in] keyword configuration keyword to find.
 * @param[in] value parsed output value.
 * @param[in] ResetLine whether to restore the original stream position after reading.
 */
template<class TValueType>
void NInputUtility::AddScalarOption
(
 std::ifstream  &inputFile,
 const char     *keyword,
 TValueType     &value,
 bool            ResetLine
)
 /*
  * Function that assigns a value from input file, otherwise exits.
  */
{
  // Initialize string.
  std::string valString;
  // Initialize line number.
  std::streampos CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringFromKeyword(inputFile, keyword, valString) )
	{
    // Output message.
		std::string message = "Keyword: ";
		message += keyword;
		message += " not found!";
		// Exit program.
		ERROR(message);
  }

  // Assign string value.
  std::istringstream istr(valString);
  istr >> value;

	// Check if any bad character is input.
	if( istr.fail() ) ERROR("Wrong character input inside expected parameter type!");

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that works as a template to read vector data, otherwise assigns default value.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in,out] inputFile input stream whose read position may change.
 * @param[in] keyword configuration keyword to find.
 * @param[in] value parsed output value.
 * @param[in] defaultValue value used when the keyword is absent.
 * @param[in] ResetLine whether to restore the original stream position after reading.
 */
template<class TValueType>
void NInputUtility::AddVectorOption
(
 std::ifstream           &inputFile,
 const char              *keyword,
 std::vector<TValueType> &value,
 std::vector<TValueType>  defaultValue,
 bool                     ResetLine
)
 /*
	* Function that searched for keyword from input file and assigns it to
	* the templated value, otherwise uses default value.
	*/
{
	// First, ensure that the value vector is empty.
	value.clear();

  // Initialize string.
	std::vector<std::string> valStringVector;
  // Initialize line number.
  std::streampos CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringVectorFromKeyword(inputFile, keyword, valStringVector) )
	{
		// Assign default value.
		value = defaultValue;
    // Go back to previous line location.
    if( ResetLine ) { inputFile.clear(); inputFile.seekg(CurrentLine); }
    return;
  }

	// Assign string value.
	for(std::size_t iData=0; iData<valStringVector.size(); iData++)
	{
		TValueType tmp;
		std::istringstream istr(valStringVector[iData]);
		istr >> tmp;

		// Check if any bad character is input.
		if( istr.fail() ) ERROR("Wrong character input inside expected parameter type!");

		// Add to input value vector.
		value.push_back(tmp);
	}

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that works as a template to read vector data, otherwise exits.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in,out] inputFile input stream whose read position may change.
 * @param[in] keyword configuration keyword to find.
 * @param[in] value parsed output value.
 * @param[in] ResetLine whether to restore the original stream position after reading.
 */
template<class TValueType>
void NInputUtility::AddVectorOption
(
 std::ifstream           &inputFile,
 const char              *keyword,
 std::vector<TValueType> &value,
 bool                     ResetLine
)
 /*
	* Function that searched for keyword from input file and assigns it to
	* the templated value, otherwise exits.
	*/
{
	// First, ensure that the value vector is empty.
	value.clear();
  
	// Initialize string.
	std::vector<std::string> valStringVector;
  // Initialize line number.
  std::streampos CurrentLine;
  // Record current line.
  if( ResetLine ) CurrentLine = inputFile.tellg();

  // Use default value, if string could not be found.
  if( !FindStringVectorFromKeyword(inputFile, keyword, valStringVector) )
	{
    // Output message.
		std::string message = "Keyword: ";
		message += keyword;
		message += " not found!";
		// Exit program.
		ERROR(message);
  }

	// Assign string value.
	for(std::size_t iData=0; iData<valStringVector.size(); iData++)
	{
		TValueType tmp;
		std::istringstream istr(valStringVector[iData]);
		istr >> tmp;

		// Check if any bad character is input.
		if( istr.fail() ) ERROR("Wrong character input inside expected parameter type!");

		// Add to input value vector.
		value.push_back(tmp);
	}

  // Go back to previous line location.
  if( ResetLine ) inputFile.seekg(CurrentLine);
}


