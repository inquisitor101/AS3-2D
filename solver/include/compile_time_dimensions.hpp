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
	APPLY( 6, 7,4) \
	APPLY( 6, 8,4) \
	APPLY( 7, 7,4) \
	APPLY( 7,10,4) \
	APPLY( 8, 8,4) \
	APPLY( 8,11,4) \
	APPLY( 9, 9,4) \
	APPLY( 9,13,4) \
	APPLY(10,10,4) \
	APPLY(10,14,4)

//-----------------------------------------------------------------------------------

/*!
 * @brief Runtime queries and compile-time limits for the factory specializations.
 *
 * @note The table is generated from AS3_SUPPORTED_ELEMENT_SPECIALIZATIONS,
 *   which also generates the factory dispatch.
 */
struct CSupportedElementDimensions
{
	/*!
	 * @brief One supported combination of solution points, integration points and variables.
	 */
	struct CEntry
	{
		std::size_t mNSol1D; ///< Number of solution points in 1D.
		std::size_t mNInt1D; ///< Number of integration points in 1D.
		std::size_t mNVar;   ///< Number of working variables.
	};

	/*!
	 * @brief All combinations instantiated by the factory.
	 */
	inline static constexpr CEntry mCombinations[] =
	{
#define AS3_DIMENSION_ENTRY(K,M,N) {K,M,N},
		AS3_SUPPORTED_ELEMENT_SPECIALIZATIONS(AS3_DIMENSION_ENTRY)
#undef AS3_DIMENSION_ENTRY
	};

	/*!
	 * @brief Getter function which returns the largest supported solution-point count in 1D.
	 *
	 * @return maximum mNSol1D in the specialization table.
	 */
	static constexpr std::size_t GetMaxSol1D(void)
	{
		std::size_t nMax = 0;
		for(const auto &entry : mCombinations)
		{
			if( entry.mNSol1D > nMax ) nMax = entry.mNSol1D;
		}
		return nMax;
	}

	/*!
	 * @brief Getter function which returns the largest supported integration-point count in 1D.
	 *
	 * @return maximum mNInt1D in the specialization table.
	 */
	static constexpr std::size_t GetMaxInt1D(void)
	{
		std::size_t nMax = 0;
		for(const auto &entry : mCombinations)
		{
			if( entry.mNInt1D > nMax ) nMax = entry.mNInt1D;
		}
		return nMax;
	}

	/*!
	 * @brief Function that checks whether a point-count pair has a factory specialization.
	 *
	 * @param[in] nSol1D number of solution points in 1D.
	 * @param[in] nInt1D number of integration points in 1D.
	 *
	 * @return true if the pair occurs with any supported variable count.
	 */
	static constexpr bool IsSupported(std::size_t nSol1D,
	                                  std::size_t nInt1D)
	{
		for(const auto &entry : mCombinations)
		{
			if( entry.mNSol1D == nSol1D && entry.mNInt1D == nInt1D ) return true;
		}
		return false;
	}

	/*!
	 * @brief Function that checks whether an exact factory specialization exists.
	 *
	 * @param[in] nSol1D number of solution points in 1D.
	 * @param[in] nInt1D number of integration points in 1D.
	 * @param[in] nVar number of working variables.
	 *
	 * @return true if the complete combination occurs in the specialization table.
	 */
	static constexpr bool IsSupported(std::size_t nSol1D,
	                                  std::size_t nInt1D,
	                                  std::size_t nVar)
	{
		for(const auto &entry : mCombinations)
		{
			if( entry.mNSol1D == nSol1D &&
			    entry.mNInt1D == nInt1D &&
			    entry.mNVar   == nVar ) return true;
		}
		return false;
	}
};
