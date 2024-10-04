#pragma once

#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <cmath>
#include "error_structure.hpp"
#include "data_structure.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif



/* * * 
 * Global name-aliasing definitions.
 * */

typedef double as3double;                                                  ///< Type definition of the floating precision.


template<typename T>
using as3vector1d = std::vector<T>;                                        ///< Vector of a data type: 1D.

template<typename T>
using as3vector2d = std::vector<std::vector<T>>;                           ///< Vector of a data type: 2D.

template<typename T>
using as3vector3d = std::vector<std::vector<std::vector<T>>>;              ///< Vector of a data type: 3D.

template<typename T>
using as3vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>; ///< Vector of a data type: 4D.




/* * * 
 * Global variable definitions.
 * */

// Abbreviation for specific AS3 convention values.
const unsigned int AS3_MAGIC_NUMBER = 3735929054; ///< AS3 file magic number.
const int          CGNS_STRING_SIZE = 33;         ///< CGNS string length.

// Abbreviation of common floating values used in the code.
const as3double C_PI   = static_cast<as3double>( 3.14159265358979323846264338327950288419 ); ///< Value of pi. 
const as3double C_ZERO = static_cast<as3double>( 0.0 );                                      ///< Value of 0.
const as3double C_HALF = static_cast<as3double>( 0.5 );                                      ///< Value of 0.5.
const as3double C_ONE  = static_cast<as3double>( 1.0 );                                      ///< Value of 1.
const as3double C_TWO  = static_cast<as3double>( 2.0 );                                      ///< Value of 2.
const as3double C_180  = static_cast<as3double>( 180.0 );                                    ///< Value of 180.

// Abbreviation for constant physical values.
const as3double C_GMA  = static_cast<as3double>(1.4);     ///< Ratio of specific heats.
const as3double C_GM1  = C_GMA - C_ONE;                   ///< Abbreviation for: gamma - 1.
const as3double C_RGAS = static_cast<as3double>(287.058); ///< Gas constant: R.
const as3double C_CV   = C_RGAS/C_GM1;                    ///< Specific heat at constant volume.
const as3double C_CP   = C_GMA*C_CV;                      ///< Specific heat at constant pressure.



/* * * 
 * Global enumeration definitions.
 * */


/*!
 * @brief Enumerated type for mesh format.
 */
enum class EMeshFormat
{
	PLOT3D,
	AS3
};

/*!
 * @brief Map for the mesh format.
 */
const std::map<std::string, EMeshFormat>
MapMeshFormat = 
{
	{ "PLOT3D", EMeshFormat::PLOT3D },
	{ "AS3",    EMeshFormat::AS3    }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the arrangement of the DOFs.
 */
enum class ETypeDOF
{
	LGL,
	EQD
};

/*!
 * @brief Map for the type of the DOFs.
 */
const std::map<std::string, ETypeDOF>
MapTypeDOF = 
{
	{ "LGL", ETypeDOF::LGL },
	{ "EQD", ETypeDOF::EQD }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the available BCs.
 */
enum class ETypeBC
{
	INTERFACE
};

/*!
 * @brief Map for the type of boundary conditions.
 */
const std::map<std::string, ETypeBC>
MapTypeBCs = 
{
	{ "INTERFACE",  ETypeBC::INTERFACE  }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the file format.
 */
enum class EFormatFile
{
	ASCII,
	BINARY
};

/*!
 * @brief Map for the file format.
 */
const std::map<std::string, EFormatFile>
MapFormatFile = 
{
	{ "ASCII",  EFormatFile::ASCII  },
	{ "BINARY", EFormatFile::BINARY }
};

//--------------------------------------

/*!
 * @brief Enumerated values for the (quadrilateral) element face indices.
 */
enum class EFaceElement 
{
	IMIN,
	IMAX,
	JMIN,
	JMAX
};

/*!
 * @brief Map for the element face values.
 */
const std::map<unsigned int, EFaceElement>
MapFaceElement = 
{
	{ 0, EFaceElement::IMIN },
	{ 1, EFaceElement::IMAX },
	{ 2, EFaceElement::JMIN },
	{ 3, EFaceElement::JMAX }
};

//--------------------------------------

/*!
 * @brief Enumerated type for output visualization format.
 */
enum class EVisualFormat
{
	VTK_LEGACY_BINARY
};

/*!
 * @brief Map for the visualization format.
 */
const std::map<std::string, EVisualFormat>
MapVisualFormat = 
{
	{ "VTK_LEGACY_BINARY", EVisualFormat::VTK_LEGACY_BINARY }
};

//--------------------------------------

/*!
 * @brief Enumerated type for temporal discretization.
 */
enum class ETemporalScheme
{
	SSPRK3
};

/*!
 * @brief Map for the temporal discretization.
 */
const std::map<std::string, ETemporalScheme>
MapTemporalScheme = 
{
	{ "SSPRK3", ETemporalScheme::SSPRK3 }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the solver.
 */
enum class ETypeSolver
{
	EE
};

/*!
 * @brief Map for the solvers.
 */
const std::map<std::string, ETypeSolver>
MapTypeSolver = 
{
	{ "EE", ETypeSolver::EE }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the initial condition.
 */
enum class ETypeIC
{
	GAUSSIAN_PRESSURE,
	ISENTROPIC_VORTEX
};

/*!
 * @brief Map for the initial conditions.
 */
const std::map<std::string, ETypeIC>
MapTypeIC = 
{
	{ "GAUSSIAN_PRESSURE", ETypeIC::GAUSSIAN_PRESSURE },
	{ "ISENTROPIC_VORTEX", ETypeIC::ISENTROPIC_VORTEX }
};

//--------------------------------------

/*!
 * @brief Enumerated type for variables to write.
 */
enum class EWriteVariable
{
	DENSITY,
	MOMENTUM,
	TOTAL_ENERGY,
	PRESSURE,
	VELOCITY,
	VORTICITY,
	MACH,
	TEMPERATURE,
	ENTROPY
};

/*!
 * @brief Map for the variables to write.
 */
const std::map<std::string, EWriteVariable>
MapWriteVariable = 
{
	{ "DENSITY",      EWriteVariable::DENSITY },
	{ "MOMENTUM",     EWriteVariable::MOMENTUM },
	{ "TOTAL_ENERGY", EWriteVariable::TOTAL_ENERGY },
	{ "PRESSURE",     EWriteVariable::PRESSURE },
	{ "VELOCITY",     EWriteVariable::VELOCITY },
	{ "VORTICITY",    EWriteVariable::VORTICITY },
	{ "MACH",         EWriteVariable::MACH },
	{ "TEMPERATURE",  EWriteVariable::TEMPERATURE },
	{ "ENTROPY",      EWriteVariable::ENTROPY }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the buffer layer.
 */
enum class ETypeBufferLayer
{
	NONE
};

/*!
 * @brief Map for the buffer layer.
 */
const std::map<std::string, ETypeBufferLayer>
MapTypeBufferLayer = 
{
	{ "NONE", ETypeBufferLayer::NONE }
};

//--------------------------------------

/*!
 * @brief Enumerated type for the Riemann solver.
 */
enum class ETypeRiemannSolver
{
	ROE
};

/*!
 * @brief Map for the Riemann solvers.
 */
const std::map<std::string, ETypeRiemannSolver>
MapTypeRiemannSolver = 
{
	{ "ROE", ETypeRiemannSolver::ROE }
};





/* * * 
 * Global function definitions.
 * */




