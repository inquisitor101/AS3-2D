#pragma once

#include "option_structure.hpp"


/*!
 * @brief A struct used for storing the parameters in boundary markers. 
 */
struct IBoundaryParamMarker
{
	std::string mName; ///< Name of the current marker.
	
	/*!
	 * @brief Virtual destructor, does nothing, but needed for polymorphism.
	 */
	virtual ~IBoundaryParamMarker() {}

	/*!
	 * @brief Pure virtual function that returns the type of BC on this marker. Must be overridden.
	 */
	virtual ETypeBC GetTypeBC(void) const = 0;
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A struct used for storing the parameters in periodic markers. 
 */
struct CPeriodicParamMarker : public IBoundaryParamMarker
{
	/*!
	 * @brief Constructor that defines the parameters of this class.
	 * 
	 * @param[in] buffer vector of strings containing the parameters.
	 */
	explicit CPeriodicParamMarker(as3vector1d<std::string> buffer);

	/*!
	 * @brief Function that returns the type of BC on this marker, which is periodic.
	 */
	ETypeBC GetTypeBC(void) const override {return ETypeBC::PERIODIC;}

	std::string              mNameMatching; ///< Name of the matching marker.
  std::array<as3double, 2> mVectorTrans;  ///< Translation vector, from I to J.
};

//-----------------------------------------------------------------------------------=

/*!
 * @brief A struct used for storing the parameters in symmetry markers. 
 */
struct CSymmetryParamMarker : public IBoundaryParamMarker
{
	/*!
	 * @brief Constructor that defines the parameters of this class.
	 * 
	 * @param[in] buffer string containing the parameters.
	 */
	explicit CSymmetryParamMarker(std::string buffer)
  {
    mName = buffer;
  }

	/*!
	 * @brief Function that returns the type of BC on this marker, which is symmetry.
	 */
	ETypeBC GetTypeBC(void) const override {return ETypeBC::SYMMETRY;}
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Parameters for a static mass flow inlet imposed through a Riemann state.
 */
struct CStaticMassFlowInletRiemannParamMarker : public IBoundaryParamMarker
{
  /*!
   * @param[in] buffer marker name, density, velocity_x and velocity_y.
   */
  explicit CStaticMassFlowInletRiemannParamMarker(const as3vector1d<std::string> &buffer);

  ETypeBC GetTypeBC(void) const override {return ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN;}

  as3double                mDensity;  ///< Prescribed static density.
  std::array<as3double, 2> mVelocity; ///< Prescribed velocity components.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Parameters for a pressure outlet imposed through a Riemann state.
 */
struct CPressureOutletRiemannParamMarker : public IBoundaryParamMarker
{
  /*!
   * @param[in] buffer marker name and pressure.
   */
  explicit CPressureOutletRiemannParamMarker(const as3vector1d<std::string> &buffer);

  ETypeBC GetTypeBC(void) const override {return ETypeBC::PRESSURE_OUTLET_RIEMANN;}

  as3double mPressure; ///< Prescribed static pressure.
};
