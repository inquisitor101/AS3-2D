#pragma once

/*!
 * @file
 * @brief Parsed parameters for physical and periodic boundary markers.
 */


#include <array>
#include <string>

#include "option_structure.hpp"


/*!
 * @brief A struct used for storing the parameters in boundary markers.
 */
struct IBoundaryParamMarker
{
	std::string mName; ///< Name of the current marker.
	

	/*!
	 * @brief Destructor that releases the resources owned by this object.
	 */
	virtual ~IBoundaryParamMarker() {}

	/*!
	 * @brief Getter function which returns the boundary-condition type implemented by this object.
	 *
	 * @return boundary-condition type implemented by this object.
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
	 * @brief Getter function which returns the boundary-condition type implemented by this object.
	 *
	 * @return boundary-condition type implemented by this object.
	 */
	ETypeBC GetTypeBC(void) const override {return ETypeBC::PERIODIC;}

	std::string              mNameMatching;  ///< Name of the matching marker.
  std::array<as3double, 2> mVectorTrans;         ///< Translation vector, from I to J.
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
	 * @brief Getter function which returns the boundary-condition type implemented by this object.
	 *
	 * @return boundary-condition type implemented by this object.
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
   * @brief Constructor that initializes CStaticMassFlowInletRiemannParamMarker from the specified data.
   *
   * @param[in] buffer marker name, density, velocity_x and velocity_y.
   */
  explicit CStaticMassFlowInletRiemannParamMarker(const as3vector1d<std::string> &buffer);

  /*!
   * @brief Getter function which returns the boundary-condition type implemented by this marker or boundary.
   *
   * @return boundary-condition type implemented by this object.
   */
  ETypeBC GetTypeBC(void) const override {return ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN;}

  as3double                mDensity;   ///< Prescribed static density.
  std::array<as3double, 2> mVelocity;  ///< Prescribed velocity components.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Parameters for a pressure outlet imposed through a Riemann state.
 */
struct CPressureOutletRiemannParamMarker : public IBoundaryParamMarker
{

  /*!
   * @brief Constructor that initializes CPressureOutletRiemannParamMarker from the specified data.
   *
   * @param[in] buffer marker name and pressure.
   */
  explicit CPressureOutletRiemannParamMarker(const as3vector1d<std::string> &buffer);

  /*!
   * @brief Getter function which returns the boundary-condition type implemented by this marker or boundary.
   *
   * @return boundary-condition type implemented by this object.
   */
  ETypeBC GetTypeBC(void) const override {return ETypeBC::PRESSURE_OUTLET_RIEMANN;}

  as3double mPressure; ///< Prescribed static pressure.
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Parameters for a farfield boundary.
 */
struct CFarfieldParamMarker : public IBoundaryParamMarker
{

  /*!
   * @brief Constructor that initializes CFarfieldParamMarker from the specified data.
   *
   * @param[in] buffer marker name, density, velocity_x, velocity_y and pressure.
   */
  explicit CFarfieldParamMarker(const as3vector1d<std::string> &buffer);

  /*!
   * @brief Getter function which returns the boundary-condition type implemented by this marker or boundary.
   *
   * @return boundary-condition type implemented by this object.
   */
  ETypeBC GetTypeBC(void) const override {return ETypeBC::FARFIELD;}

  as3double                mDensity;   ///< Prescribed static density.
  std::array<as3double, 2> mVelocity;  ///< Prescribed velocity components.
  as3double                mPressure;  ///< Prescribed static pressure.
};
