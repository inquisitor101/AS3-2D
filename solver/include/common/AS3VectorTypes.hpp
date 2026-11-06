#pragma once

/*!
 * @file
 * @brief Solver aliases for standard vector containers.
 */


#include <vector>


/*!
 * @brief Alias for as3vector1d.
 */
template<typename T>
using as3vector1d = std::vector<T>;                                        ///< Vector of a data type: 1D.

/*!
 * @brief Alias for as3vector2d.
 */
template<typename T>
using as3vector2d = std::vector<std::vector<T>>;                           ///< Vector of a data type: 2D.

/*!
 * @brief Alias for as3vector3d.
 */
template<typename T>
using as3vector3d = std::vector<std::vector<std::vector<T>>>;              ///< Vector of a data type: 3D.

/*!
 * @brief Alias for as3vector4d.
 */
template<typename T>
using as3vector4d = std::vector<std::vector<std::vector<std::vector<T>>>>; ///< Vector of a data type: 4D.
