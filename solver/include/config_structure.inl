

//-----------------------------------------------------------------------------------
// Implementation of the templated functions in CConfig. 
//-----------------------------------------------------------------------------------

template<class TValueType>
void CConfig::PadEntriesVectorData
(
 std::vector<TValueType> &data,
 std::string        			keyword,
 unsigned short           nExpected,
 unsigned short           nCondition1,
 unsigned short           nCondition2,
 unsigned short           nCondition3
)
 /*
  * Function that pads the remaining entries in a given vector of templated value.
  */
{
  // Size of current data vector.
  size_t nData = data.size();

  // Check first condition.
  if( nData != nExpected )
	{
    // Data flag.
    bool PadData = false;
    // Check if there need be any modification, if so flag the data.
    if( nData == nCondition1 ) PadData = true;
    if( nData == nCondition2 ) PadData = true;
    if( nData == nCondition3 ) PadData = true;

    // Pad data accordingly.
    if( PadData ) 
			for(size_t i=nData; i<nExpected; i++) data.push_back( data[i-1] );
  }

  // Consistency check.
  if( data.size() != nExpected )
	{
    std::string message = "Keyword: ";
    message += keyword;
    message += " must be of size: ";
    message += std::to_string(nExpected);
    ERROR(message);
  }
}

//-----------------------------------------------------------------------------------

template<class TValueType>
void CConfig::PadEntriesVectorDefaultData
(
 std::vector<TValueType> &data,
 std::vector<TValueType> &reference,
 std::string        			keyword
)
 /*
  * Function that pads the remaining entries in a given vector of templated value
  * based on a default reference value.
  */
{
  // Size of current data vector.
  unsigned short nData     = data.size();
  // Size of expected default reference data vector.
  unsigned short nExpected = reference.size();

  // Check if padding is needed.
  if( nData != nExpected ) 
		for(unsigned short i=nData; i<nExpected; i++) data.push_back( reference[i] );

  // Consistency check.
  if( data.size() != nExpected )
	{
    std::string message = "Keyword: ";
    message += keyword;
    message += " must be of size: ";
    message += std::to_string(nExpected);
    ERROR(message);
  }
}

//-----------------------------------------------------------------------------------

template<class Tenum>
Tenum CConfig::GenericScalarMap
(
 const std::map<std::string, Tenum> &mapper,
 std::string                        &key,
 std::string                         msg
)
 /*
	* Function that maps a single string to an enum.
	*/
{
	// Check if data abides by map convention.
	if( !mapper.contains(key) ) ERROR(msg + " could not be mapped.");

	// Assign data according to dedicated enum.
	return mapper.at( key );
}

//-----------------------------------------------------------------------------------

template<class Tenum>
as3vector1d<Tenum> CConfig::GenericVectorMap
(
 const std::map<std::string, Tenum> &mapper,
 as3vector1d<std::string>           &key,
 std::string                         msg
)
 /*
	* Function that maps a vector of strings to an enum.
	*/
{
	as3vector1d<Tenum> val(key.size());
	for(size_t i=0; i<key.size(); i++)
	{
		val[i] = GenericScalarMap(mapper, key[i], msg);
	}
	return val;
}
