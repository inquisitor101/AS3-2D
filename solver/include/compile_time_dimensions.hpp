#pragma once

/*!
 * @file
 * @brief Compile-time dimension descriptors and supported factory specializations.
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

//-----------------------------------------------------------------------------------

/*!
 * @brief Apply a macro to each supported factory combination of (K,M,N).
 *
 * @param[in] APPLY macro accepting the solution-point count, integration-point count and variable count.
 * @note This list controls factory dispatch; the dimension descriptor itself can represent other combinations.
 */
#define AS3_SUPPORTED_ELEMENT_SPECIALIZATIONS(APPLY) \
	APPLY( 2, 2,4) \
	APPLY( 2, 3,4) \
	APPLY( 3, 3,4) \
	APPLY( 3, 4,4) \
	APPLY( 4, 4,4) \
	APPLY( 4, 5,4) \
	APPLY( 5, 5,4) \
	APPLY( 5, 7,4) \
	APPLY( 6, 6,4) \
	APPLY( 6, 8,4) \
	APPLY( 7, 7,4) \
	APPLY( 7,10,4) \
	APPLY( 8, 8,4) \
	APPLY( 8,11,4) \
	APPLY( 9, 9,4) \
	APPLY( 9,13,4) \
	APPLY(10,10,4) \
	APPLY(10,14,4)
