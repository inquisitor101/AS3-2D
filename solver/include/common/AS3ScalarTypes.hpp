#pragma once

/*!
 * @file
 * @brief Solver floating-point precision and matrix scalar constraint.
 */


#include <concepts>
#include <type_traits>


typedef double as3double; ///< Floating-point precision used by the solver.


/*!
 * @brief Mutable built-in arithmetic scalars supported by all AS3 matrix containers.
 *
 * These types have numeric zero value-initialization, non-throwing element lifetimes,
 * and bytewise copy semantics. Trivially copyable alone would also admit pointers
 * and custom classes with throwing constructors. Const/volatile entries cannot
 * satisfy the containers' ordinary mutable storage and memcpy contract.
 * This constraint is checked at compile time and adds no runtime work.
 */
template<typename T>
concept AS3MatrixScalar = (std::integral<T> || std::floating_point<T>) &&
                           std::same_as<T, std::remove_cv_t<T>>;
