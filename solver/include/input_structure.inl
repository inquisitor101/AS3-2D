

//-----------------------------------------------------------------------------------
// Implementation of the templated functions in NInputUtility.
//-----------------------------------------------------------------------------------


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
	for(unsigned short i=0; i<str.size(); i++)
	{
		val[i] = ExtractLastDigitFromString<TValueType>(str[i], sym);
	}
}

//-----------------------------------------------------------------------------------

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
	size_t pos = str.find(sym);
	
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
  size_t CurrentLine;
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
  size_t CurrentLine;
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
  size_t CurrentLine;
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
	for(unsigned short iData=0; iData<valStringVector.size(); iData++)
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
  size_t CurrentLine;
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
	for(unsigned short iData=0; iData<valStringVector.size(); iData++)
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


