#pragma once

/*!
 * @file
 * @brief Compile-time dimensions for tensor operations.
 */


#include <cstddef>


/*!
 * @brief Compile-time dimensions for a given number of solution points, integration points and variables.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 */
template<std::size_t K, std::size_t M, std::size_t N>
struct CCompileTimeDimensions
{
	inline static constexpr std::size_t nSol1D = K;    ///< Number of solution points in 1D.
	inline static constexpr std::size_t nInt1D = M;    ///< Number of integration points in 1D.
	inline static constexpr std::size_t nSol2D = K*K;  ///< Number of solution points in 2D.
	inline static constexpr std::size_t nInt2D = M*M;  ///< Number of integration points in 2D.
	inline static constexpr std::size_t nVar   = N;    ///< Number of working variables.
};
