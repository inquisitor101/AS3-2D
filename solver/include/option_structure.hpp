#pragma once

#include <sstream>
#include <vector>
#include <memory>
#include <map>
#include <cmath>
#include "error_structure.hpp"
#include "data_structure.hpp"
#include "log_structure.hpp"


/* * * 
 * Global name-aliasing definitions.
 * */

typedef double as3double;                                                  ///< Type definition of double precision.


template<typename T>
using as3vector1d = std::vector<T>;                                        ///< Vector of data type: 1D.

template<typename T>
using as3vector2d = std::vector<std::vector<T>>;                           ///< Vector of data type: 2D.

template<typename T>
using as3vector3d = std::vector<std::vector<std::vector<T>>>;              ///< Vector of data type: 3D.

template<typename T>
using as3vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>; ///< Vector of data type: 4D.




/* * * 
 * Global variable definitions.
 * */

const unsigned int AS3_MAGIC_NUMBER = 3735929054;                   ///< AS3 file magic number.
const int          CGNS_STRING_SIZE = 33;                           ///< CGNS string length.
const as3double PI_CONSTANT         = 4.0*std::atan(1.0);                ///< Value of the constant \pi.
const as3double GAMMA               = 1.4;                          ///< Ratio of specific heats.
const as3double GAMMA_MINUS_ONE     = GAMMA-1.0;                    ///< Abbreviation for: gamma - 1.
const as3double GAS_CONSTANT        = 287.058;                      ///< Gas constant: R.
const as3double CV_CONSTANT         = GAS_CONSTANT/GAMMA_MINUS_ONE; ///< Specific heat at constant volume.
const as3double CP_CONSTANT         = GAMMA*CV_CONSTANT;            ///< Specific heat at constant pressure.




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
enum class ETypeBCs
{
	PERIODIC,
	INTERFACE
};

/*!
 * @brief Map for the type of boundary conditions.
 */
const std::map<std::string, ETypeBCs>
MapTypeBCs = 
{
	{ "PERIODIC",  ETypeBCs::PERIODIC  },
	{ "INTERFACE", ETypeBCs::INTERFACE }
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
 * @brief Enumerated type for type of solver.
 */
enum class ETypeSolver
{
	EE
};

/*!
 * @brief Map for the type of solver.
 */
const std::map<std::string, ETypeSolver>
MapTypeSolver = 
{
	{ "EE", ETypeSolver::EE }
};

//--------------------------------------

/*!
 * @brief Enumerated type for type of IC.
 */
enum class ETypeIC
{
	GAUSSIAN_PRESSURE,
	ISENTROPIC_VORTEX
};

/*!
 * @brief Map for the type of IC.
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
 * @brief Enumerated type for type of buffer layer.
 */
enum class ETypeBufferLayer
{
	NONE
};

/*!
 * @brief Map for the type of buffer layer.
 */
const std::map<std::string, ETypeBufferLayer>
MapTypeBufferLayer = 
{
	{ "NONE", ETypeBufferLayer::NONE }
};





/* * * 
 * Global function definitions.
 * */




