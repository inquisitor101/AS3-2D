#pragma once

/*!
 * @file
 * @brief Shared numerical, physical and file-format constants.
 */


#include <algorithm>
#include <cstdint>
#include <limits>

#include "common/AS3ScalarTypes.hpp"


/* * * 
 * Global variable definitions.
 * */

// Abbreviation for specific AS3 convention values.
inline constexpr std::uint32_t AS3_MAGIC_NUMBER = 3735929054;  ///< AS3 file magic number.
inline constexpr int          CGNS_STRING_SIZE = 33;           ///< CGNS string length.

// Abbreviation of common floating values used in the code.
inline constexpr as3double C_PI    = static_cast<as3double>( 3.14159265358979323846 );  ///< Value of pi. 
inline constexpr as3double C_ZERO  = static_cast<as3double>( 0.0 );                     ///< Value of 0.
inline constexpr as3double C_HALF  = static_cast<as3double>( 0.5 );                     ///< Value of 0.5.
inline constexpr as3double C_ONE   = static_cast<as3double>( 1.0 );                     ///< Value of 1.
inline constexpr as3double C_TWO   = static_cast<as3double>( 2.0 );                     ///< Value of 2.
inline constexpr as3double C_180   = static_cast<as3double>( 180.0 );                   ///< Value of 180.

// Precision and absolute tolerance for dimensionless reference-element checks.
inline constexpr as3double C_EPS            = std::numeric_limits<as3double>::epsilon();         ///< Floating-point spacing just above one.
inline constexpr as3double C_REFERENCE_ATOL = std::max(static_cast<as3double>(1.0e-10), C_EPS);  ///< Absolute floor for reference-element checks.


// Abbreviation for constant physical values.
inline constexpr as3double C_GMA   = static_cast<as3double>(1.4);      ///< Ratio of specific heats.
inline constexpr as3double C_GM1   = C_GMA - C_ONE;                    ///< Abbreviation for: gamma - 1.
inline constexpr as3double C_OVGM1 = C_ONE / C_GM1;                    ///< Abbreviation for: 1 / (gamma-1)
inline constexpr as3double C_OVGMA = C_ONE / C_GMA;                    ///< Abbreviation for: 1 / gamma
inline constexpr as3double C_RGAS  = static_cast<as3double>(287.058);  ///< Gas constant: R.
inline constexpr as3double C_CV    = C_RGAS/C_GM1;                     ///< Specific heat at constant volume.
inline constexpr as3double C_CP    = C_GMA*C_CV;                       ///< Specific heat at constant pressure.
