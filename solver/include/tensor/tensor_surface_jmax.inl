/*!
 * @file
 * @brief Compile-time tensor kernels for surface jmax.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the JMAX surface templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the JMAX surface
 *   integration points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B solution values [N][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::CompileTimeSurfaceJMAX
(
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that computes the tensor product on the surface integration points of JMAX.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A        = mLagrangeInt1DTrans.data();
	const as3double *ADer     = mDerLagrangeInt1DTrans.data();
	const as3double *aDerFace = mDerLagrangeMaxFace1D.data();

	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) solution points in 2D and abbreviation.
	constexpr std::size_t K2   = K*K;
	constexpr std::size_t K2mK = K2-K;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
		// Create temporary storage.
		as3double tmpI[M];
		as3double tmpJ[K];

		// Store boundary surface info in tmpJ.
    for(std::size_t s=0; s<K; s++) tmpJ[s] = B[l*K2+s+K2mK];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*M];

			for(std::size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t k=0; k<M; k++)
			{
        for(std::size_t ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to C.
			for(std::size_t i=0; i<M; i++) c[i] = tmpI[i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*M];

			for(std::size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t k=0; k<M; k++)
			{
        for(std::size_t ii=0; ii<K; ii++)
					tmpI[k] += aDer[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerR.
			for(std::size_t i=0; i<M; i++) cDerR[i] = tmpI[i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*M];

			// Cast array B from 1D to 2D.
			const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];
			
			for(std::size_t s=0; s<K; s++) tmpJ[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t i=0; i<K; i++)
			{
        for(std::size_t jj=0; jj<K; jj++)
					tmpJ[i] += aDerFace[jj] * b[jj][i];
			}

			for(std::size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t k=0; k<M; k++)
			{
        for(std::size_t ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerS.
			for(std::size_t i=0; i<M; i++) cDerS[i] = tmpI[i];
		}
	}
}


//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the JMAX surface
 *   integration points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] n number of variables processed with the fixed K/M operators.
 * @param[in] B solution values [n][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::PartialCompileTimeSurfaceJMAX
(
 std::size_t           n,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that computes the tensor product on the surface integration points of JMAX.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A        = mLagrangeInt1DTrans.data();
	const as3double *ADer     = mDerLagrangeInt1DTrans.data();
	const as3double *aDerFace = mDerLagrangeMaxFace1D.data();

	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) solution points in 2D and abbreviation.
	constexpr std::size_t K2   = K*K;
	constexpr std::size_t K2mK = K2-K;


	// Loop over each n entry.
	for(std::size_t l=0; l<n; l++)
	{
		// Create temporary storage.
		as3double tmpI[M];
		as3double tmpJ[K];

		// Store boundary surface info in tmpJ.
    for(std::size_t s=0; s<K; s++) tmpJ[s] = B[l*K2+s+K2mK];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*M];

			for(std::size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t k=0; k<M; k++)
			{
        for(std::size_t ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to C.
			for(std::size_t i=0; i<M; i++) c[i] = tmpI[i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*M];

			for(std::size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t k=0; k<M; k++)
			{
        for(std::size_t ii=0; ii<K; ii++)
					tmpI[k] += aDer[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerR.
			for(std::size_t i=0; i<M; i++) cDerR[i] = tmpI[i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*M];

			// Cast array B from 1D to 2D.
			const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];
			
			for(std::size_t s=0; s<K; s++) tmpJ[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t i=0; i<K; i++)
			{
        for(std::size_t jj=0; jj<K; jj++)
					tmpJ[i] += aDerFace[jj] * b[jj][i];
			}

			for(std::size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
#pragma omp simd
			for(std::size_t k=0; k<M; k++)
			{
        for(std::size_t ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerS.
			for(std::size_t i=0; i<M; i++) cDerS[i] = tmpI[i];
		}
	}
}
