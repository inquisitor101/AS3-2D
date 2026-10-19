#include "param_structure.hpp" 


//-----------------------------------------------------------------------------------
// CPeriodicParamMarker member functions.
//-----------------------------------------------------------------------------------
namespace
{
  // Convert one complete token, rejecting invalid numbers and trailing characters.
  as3double ParseBoundaryValue(const std::string &token,
                               const std::string &markerName,
                               const char        *field)
  {
    as3double value;
    std::istringstream input(token);
    if( !(input >> value) || !std::isfinite(value) )
      ERROR("Invalid " + std::string(field) + " for boundary marker: " + markerName);
  
    input >> std::ws;
    if( !input.eof() )
      ERROR("Invalid " + std::string(field) + " for boundary marker: " + markerName);
  
    return value;
  }
} // namespace


//-----------------------------------------------------------------------------------
// CPeriodicParamMarker member functions.
//-----------------------------------------------------------------------------------


CPeriodicParamMarker::CPeriodicParamMarker
(
 as3vector1d<std::string> buffer
)
 /*
	* Constructor that initializes the parameters needed in an periodic BC.
	*/
{
	// Deduce the translation vector from the buffer string and store it into tmp.
	as3double tmp[2];
	std::istringstream istr( buffer[2] + " " + buffer[3] );
	istr >> tmp[0] >> tmp[1];

	// Check if any bad character is input.
	if( istr.fail() ) ERROR("Wrong character input inside translation vector.");

	// This marker is from i's perspective: i to j.
	mName           = buffer[0];
	mNameMatching   = buffer[1];
	mVectorTrans[0] = tmp[0];
	mVectorTrans[1] = tmp[1];

	// Ensure the boundary markers are unique.
	if( mName == mNameMatching ) ERROR("Periodic markers must be unique.");
}

//-----------------------------------------------------------------------------------

CStaticMassFlowInletRiemannParamMarker::CStaticMassFlowInletRiemannParamMarker
(
 const as3vector1d<std::string> &buffer
)
{
  if( buffer.size() != 4 )
    ERROR("MARKER_BC_STATIC_MASS_FLOW_INLET_RIEMANN must specify (name, density, velocity_x, velocity_y).");

  mName        = buffer[0];
  mDensity     = ParseBoundaryValue(buffer[1], mName, "density");
  mVelocity[0] = ParseBoundaryValue(buffer[2], mName, "velocity_x");
  mVelocity[1] = ParseBoundaryValue(buffer[3], mName, "velocity_y");
}

//-----------------------------------------------------------------------------------

CPressureOutletRiemannParamMarker::CPressureOutletRiemannParamMarker
(
 const as3vector1d<std::string> &buffer
)
{
  if( buffer.size() != 2 )
    ERROR("MARKER_BC_PRESSURE_OUTLET_RIEMANN must specify (name, pressure).");

  mName     = buffer[0];
  mPressure = ParseBoundaryValue(buffer[1], mName, "pressure");
}

