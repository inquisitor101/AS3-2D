/*!
 * @file
 * @brief Inline implementations for simulation configuration, parameter lookup and option validation.
 */

#include <cstddef>
#include <map>
#include <string>
#include <vector>



//-----------------------------------------------------------------------------------
// Implementation of the templated functions in CConfig. 
//-----------------------------------------------------------------------------------

/*!
 * @brief Function that pads remaining entries in a vector.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in,out] data data vector to pad.
 * @param[in] keyword configuration keyword to find.
 * @param[in] nExpected required number of entries.
 * @param[in] nCondition1 an allowed input length that can be padded.
 * @param[in] nCondition2 an allowed input length that can be padded.
 * @param[in] nCondition3 an allowed input length that can be padded.
 */
template<class TValueType>
void CConfig::PadEntriesVectorData
(
 std::vector<TValueType> &data,
 std::string        			keyword,
 std::size_t           nExpected,
 std::size_t           nCondition1,
 std::size_t           nCondition2,
 std::size_t           nCondition3
)
 /*
  * Function that pads the remaining entries in a given vector of templated value.
  */
{
  // Size of current data vector.
  std::size_t nData = data.size();

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
			for(std::size_t i=nData; i<nExpected; i++) data.push_back( data[i-1] );
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

/*!
 * @brief Function that pads remaining entries in a vector based on default reference data.
 *
 * @tparam TValueType type used for parsed or padded values.
 *
 * @param[in,out] data data vector to pad.
 * @param[in] reference default entries used to pad the data.
 * @param[in] keyword configuration keyword to find.
 */
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
  std::size_t nData     = data.size();
  // Size of expected default reference data vector.
  std::size_t nExpected = reference.size();

  // Check if padding is needed.
  if( nData != nExpected ) 
		for(std::size_t i=nData; i<nExpected; i++) data.push_back( reference[i] );

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

/*!
 * @brief Function that maps a single string to an enum.
 *
 * @tparam Tenum enumeration type associated with the keyword map.
 *
 * @param[in] mapper mapping from keyword strings to enum values.
 * @param[in] key keyword or keywords to map.
 * @param[in] msg prefix used in an error message.
 */
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

/*!
 * @brief Function that maps a vector of strings to an enum.
 *
 * @tparam Tenum enumeration type associated with the keyword map.
 *
 * @param[in] mapper mapping from keyword strings to enum values.
 * @param[in] key keyword or keywords to map.
 * @param[in] msg prefix used in an error message.
 */
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
	for(std::size_t i=0; i<key.size(); i++)
	{
		val[i] = GenericScalarMap(mapper, key[i], msg);
	}
	return val;
}
