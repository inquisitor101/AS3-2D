#pragma once

/*!
 * @file
 * @brief Solver scalar types, constants, enumerations and keyword maps.
 */


#include <cstddef>
#include <cstdint>
#include <string>

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



/*!
 * @brief * Global name-aliasing definitions.
 */

typedef double as3double;                                                  ///< Type definition of the floating precision.


/*!
 * @brief Alias for .
 */
template<typename T>

/*!
 * @brief Alias for as3vector1d.
 */
using as3vector1d = std::vector<T>;                                        ///< Vector of a data type: 1D.

/*!
 * @brief Alias for .
 */
template<typename T>

/*!
 * @brief Alias for as3vector2d.
 */
using as3vector2d = std::vector<std::vector<T>>;                           ///< Vector of a data type: 2D.

/*!
 * @brief Alias for .
 */
template<typename T>

/*!
 * @brief Alias for as3vector3d.
 */
using as3vector3d = std::vector<std::vector<std::vector<T>>>;              ///< Vector of a data type: 3D.

/*!
 * @brief Alias for .
 */
template<typename T>

/*!
 * @brief Alias for as3vector4d.
 */
using as3vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>; ///< Vector of a data type: 4D.




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


// Abbreviation for constant physical values.
inline constexpr as3double C_GMA   = static_cast<as3double>(1.4);      ///< Ratio of specific heats.
inline constexpr as3double C_GM1   = C_GMA - C_ONE;                    ///< Abbreviation for: gamma - 1.
inline constexpr as3double C_OVGM1 = C_ONE / C_GM1;                    ///< Abbreviation for: 1 / (gamma-1)
inline constexpr as3double C_OVGMA = C_ONE / C_GMA;                    ///< Abbreviation for: 1 / gamma
inline constexpr as3double C_RGAS  = static_cast<as3double>(287.058);  ///< Gas constant: R.
inline constexpr as3double C_CV    = C_RGAS/C_GM1;                     ///< Specific heat at constant volume.
inline constexpr as3double C_CP    = C_GMA*C_CV;                       ///< Specific heat at constant pressure.



/* * * 
 * Global enumeration definitions.
 * */


/*!
 * @brief Enumerated type for mesh format.
 */
enum class EMeshFormat
{
	PLOT3D,  ///< Plot3D mesh format.
	AS3      ///< Native AS3 binary mesh format.
};

/*!
 * @brief Map for the mesh format.
 */
inline const std::map<std::string, EMeshFormat>
MapMeshFormat = 
{
	{ "PLOT3D", EMeshFormat::PLOT3D },
	{ "AS3",    EMeshFormat::AS3    }
}; ///< Map for the mesh format.

//--------------------------------------

/*!
 * @brief Enumerated type for the arrangement of the DOFs.
 */
enum class ETypeDOF
{
	LGL,  ///< Legendre-Gauss-Lobatto nodes.
	EQD   ///< Equidistant nodes.
};

/*!
 * @brief Map for the type of the DOFs.
 */
inline const std::map<std::string, ETypeDOF>
MapTypeDOF = 
{
	{ "LGL", ETypeDOF::LGL },
	{ "EQD", ETypeDOF::EQD }
}; ///< Map for the type of the DOFs.

//--------------------------------------

/*!
 * @brief Enumerated type for the supported quadratures.
 */
enum class ETypeQuadrature
{
	LEGENDRE_GAUSS_LOBATTO,  ///< Legendre-Gauss-Lobatto.
	GAUSS_LEGENDRE           ///< Gauss-Legendre.
};

/*!
 * @brief Map for the type of quadrature.
 */
inline const std::map<std::string, ETypeQuadrature>
MapTypeQuadrature = 
{
	{ "LEGENDRE_GAUSS_LOBATTO", ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO },
	{ "GAUSS_LEGENDRE",         ETypeQuadrature::GAUSS_LEGENDRE }
}; ///< Map for the type of the quadratures.

//--------------------------------------

/*!
 * @brief Enumerated type for the supported DG methods.
 */
enum class ETypeDGMethod
{
	NODAL_DG, ///< Standard nodal DG.
  DG_SEM    ///< Collocated: Legendre-Gauss-Lobatto.
};

/*!
 * @brief Map for the type of DG method.
 */
inline const std::map<std::string, ETypeDGMethod>
MapTypeDGMethod = 
{
	{ "NODAL_DG", ETypeDGMethod::NODAL_DG },
	{ "DG_SEM",   ETypeDGMethod::DG_SEM }
}; ///< Map for the type of the DG methods.

//--------------------------------------

/*!
 * @brief Enumerated type for the available BCs.
 */
enum class ETypeBC
{
  PERIODIC,                        ///< Periodically paired zone edges.
  SYMMETRY,                        ///< Inviscid symmetry boundary.
  FARFIELD,                        ///< Farfield boundary.
  STATIC_MASS_FLOW_INLET_RIEMANN,  ///< Prescribed density and velocity supplied as an inlet Riemann state.
  PRESSURE_OUTLET_RIEMANN          ///< Prescribed static pressure supplied as an outlet Riemann state.
};

/*!
 * @brief Map for the type of boundary conditions.
 */
inline const std::map<std::string, ETypeBC>
MapTypeBCs = 
{
	{ "PERIODIC",  ETypeBC::PERIODIC },
  { "SYMMETRY", ETypeBC::SYMMETRY },
  { "FARFIELD", ETypeBC::FARFIELD },
  { "STATIC_MASS_FLOW_INLET_RIEMANN", ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN },
  { "PRESSURE_OUTLET_RIEMANN", ETypeBC::PRESSURE_OUTLET_RIEMANN }
}; ///< Map for the type of boundary conditions.

//--------------------------------------

/*!
 * @brief Enumerated type for the file format.
 */
enum class EFormatFile
{
	ASCII,  ///< Text representation.
	BINARY  ///< Binary representation.
};

/*!
 * @brief Map for the file format.
 */
inline const std::map<std::string, EFormatFile>
MapFormatFile = 
{
	{ "ASCII",  EFormatFile::ASCII  },
	{ "BINARY", EFormatFile::BINARY }
}; ///< Map for the file format.

//--------------------------------------

/*!
 * @brief Enumerated values for the (quadrilateral) element face indices.
 */
enum class EFaceLocation 
{
	IMIN,  ///< Minimum local i-edge.
	IMAX,  ///< Maximum local i-edge.
	JMIN,  ///< Minimum local j-edge.
	JMAX   ///< Maximum local j-edge.
};

/*!
 * @brief Map for the element face values.
 */
inline const std::map<std::size_t, EFaceLocation>
MapFaceElement = 
{
	{ 0, EFaceLocation::IMIN },
	{ 1, EFaceLocation::IMAX },
	{ 2, EFaceLocation::JMIN },
	{ 3, EFaceLocation::JMAX }
}; ///< Map for the element face values.

//--------------------------------------

/*!
 * @brief Enumerated type for output visualization format.
 */
enum class EVisualFormat
{
	VTK_LEGACY_BINARY ///< Legacy VTK file with binary payload.
};

/*!
 * @brief Map for the visualization format.
 */
inline const std::map<std::string, EVisualFormat>
MapVisualFormat = 
{
	{ "VTK_LEGACY_BINARY", EVisualFormat::VTK_LEGACY_BINARY }
}; ///< Map for the visualization format.

//--------------------------------------

/*!
 * @brief Enumerated type for temporal discretization.
 */
enum class ETemporalScheme
{
	SSPRK3 ///< Three-stage, third-order strong-stability-preserving Runge-Kutta method.
};

/*!
 * @brief Map for the temporal discretization.
 */
inline const std::map<std::string, ETemporalScheme>
MapTemporalScheme = 
{
	{ "SSPRK3", ETemporalScheme::SSPRK3 }
}; ///< Map for the temporal discretization.

//--------------------------------------

/*!
 * @brief Enumerated type for the solver.
 */
enum class ETypeSolver
{
	EE ///< Compressible Euler equations.
};

/*!
 * @brief Map for the solvers.
 */
inline const std::map<std::string, ETypeSolver>
MapTypeSolver = 
{
	{ "EE", ETypeSolver::EE }
}; ///< Map for the solvers.

//--------------------------------------

/*!
 * @brief Enumerated type for the initial condition.
 */
enum class ETypeIC
{
	GAUSSIAN_PRESSURE,  ///< Gaussian pressure disturbance.
	ISENTROPIC_VORTEX   ///< Isentropic vortex disturbance.
};

/*!
 * @brief Map for the initial conditions.
 */
inline const std::map<std::string, ETypeIC>
MapTypeIC = 
{
	{ "GAUSSIAN_PRESSURE", ETypeIC::GAUSSIAN_PRESSURE },
	{ "ISENTROPIC_VORTEX", ETypeIC::ISENTROPIC_VORTEX }
}; ///< Map for the initial conditions.

//--------------------------------------

/*!
 * @brief Enumerated type for variables to write.
 */
enum class EWriteVariable
{
	DENSITY,       ///< Density.
	MOMENTUM,      ///< Cartesian momentum components.
	TOTAL_ENERGY,  ///< Total energy per unit volume.
	PRESSURE,      ///< Static pressure.
	VELOCITY,      ///< Cartesian velocity components.
	VORTICITY,     ///< Out-of-plane vorticity.
	MACH,          ///< Mach number.
	TEMPERATURE,   ///< Static temperature.
	ENTROPY        ///< Specific entropy output.
};

/*!
 * @brief Map for the variables to write.
 */
inline const std::map<std::string, EWriteVariable>
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
}; ///< Map for the variables to write.

//--------------------------------------

/*!
 * @brief Enumerated type for the buffer layer.
 */
enum class ETypeBufferLayer
{
	NONE ///< No buffer layer.
};

/*!
 * @brief Map for the buffer layer.
 */
inline const std::map<std::string, ETypeBufferLayer>
MapTypeBufferLayer = 
{
	{ "NONE", ETypeBufferLayer::NONE }
}; ///< Map for the buffer layer.

//--------------------------------------

/*!
 * @brief Enumerated type for the Riemann solver.
 */
enum class ETypeRiemannSolver
{
	ROE ///< Roe approximate Riemann solver.
};

/*!
 * @brief Map for the Riemann solvers.
 */
inline const std::map<std::string, ETypeRiemannSolver>
MapTypeRiemannSolver = 
{
	{ "ROE", ETypeRiemannSolver::ROE }
}; ///< Map for the Riemann solvers.


//--------------------------------------

/*!
 * @brief Enumerated type for the type of the element faces.
 */
enum class ETypeFaceGeometry
{
	INTERNAL,  ///< Face shared by two elements in the same zone.
	BOUNDARY,  ///< Physical boundary face.
  INTERFACE        ///< Face connecting two zones, including periodic pairs.
};

//--------------------------------------

/*!
 * @brief Enumerated type for the face type.
 */
enum class ETypeFace
{
  IFACE,  ///< A face at constant local i.
  JFACE   ///< A face at constant local j.
};

//--------------------------------------

/*!
 * @brief Strategy used to construct the face permutation.
 */
enum class EFaceLoadBalanceStrategy
{
  GREEDY ///< Greedy distribution using relative face costs.
};

//--------------------------------------

/*!
 * @brief Strategy used to construct the type of zone edge markers.
 */
enum class ETypeZoneMarker
{
  INTERNAL,  ///< Internal zone-edge marker connecting two zones.
  EXTERNAL   ///< Physical zone-edge marker.
};


/* * * 
 * Global function definitions.
 * */




