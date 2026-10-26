/*!
 * @file
 * @brief Compile-time tensor kernels for volume.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the volume templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the volume integration
 *   points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B solution values [N][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [N][M*M]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [N][M*M]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [N][M*M]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::CompileTimeVolume
(
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that computes the tensor product on the volume integration points.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A    = mLagrangeInt1DTrans.data();
	const as3double *ADer = mDerLagrangeInt1DTrans.data();

	// Cast the arrays from 1D to 2D, for convenience.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr std::size_t M2 = M*M; 
	constexpr std::size_t K2 = K*K;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
		// Create temporary storage.
		as3double tmpJ[K][M];
		as3double tmpI[M][M];

		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];

		for(std::size_t i=0; i<K; i++)
		{
			for(std::size_t s=0; s<M; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
			for(std::size_t jj=0; jj<K; jj++)
			{
#pragma omp simd
        for(std::size_t j=0; j<M; j++)
					tmpJ[i][j] += a[jj][j] * b[jj][i];
			}
		}
	

		// Check if the interpolated solution is needed.
		if( C )
		{
			// Cast array C from 1D to 2D.
			as3double (*c)[M] = (as3double (*)[M]) &C[l*M2];

			for(std::size_t j=0; j<M; j++)
			{
        for(std::size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
          for(std::size_t i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to C.
			for(std::size_t j=0; j<M; j++)
				for(std::size_t i=0; i<M; i++)
					c[j][i] = tmpI[j][i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Cast array CDerR from 1D to 2D.
			as3double (*cDerR)[M] = (as3double (*)[M]) &CDerR[l*M2];

			for(std::size_t j=0; j<M; j++)
			{
				for(std::size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
          for(std::size_t i=0; i<M; i++)
						tmpI[j][i] += aDer[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerR.
			for(std::size_t j=0; j<M; j++)
				for(std::size_t i=0; i<M; i++)
					cDerR[j][i] = tmpI[j][i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Cast array CDerS from 1D to 2D. 
			as3double (*cDerS)[M] = (as3double (*)[M]) &CDerS[l*M2];

			for(std::size_t i=0; i<K; i++)
			{
        for(std::size_t s=0; s<M; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<K; jj++)
				{
#pragma omp simd
          for(std::size_t j=0; j<M; j++)
						tmpJ[i][j] += aDer[jj][j] * b[jj][i];
				}
			}

			for(std::size_t j=0; j<M; j++)
			{
        for(std::size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
          for(std::size_t i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerS.
			for(std::size_t j=0; j<M; j++)
				for(std::size_t i=0; i<M; i++)
					cDerS[j][i] = tmpI[j][i];
		}
	}
}


//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the volume integration
 *   points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] n number of variables processed with the fixed K/M operators.
 * @param[in] B solution values [n][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [n][M*M]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [n][M*M]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [n][M*M]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::PartialCompileTimeVolume
(
 std::size_t           n,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that computes the tensor product on the volume integration points.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A    = mLagrangeInt1DTrans.data();
	const as3double *ADer = mDerLagrangeInt1DTrans.data();

	// Cast the arrays from 1D to 2D, for convenience.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr std::size_t M2 = M*M; 
	constexpr std::size_t K2 = K*K;


	// Loop over each n entry.
	for(std::size_t l=0; l<n; l++)
	{
		// Create temporary storage.
		as3double tmpJ[K][M];
		as3double tmpI[M][M];

		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];

		for(std::size_t i=0; i<K; i++)
		{
			for(std::size_t s=0; s<M; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
			for(std::size_t jj=0; jj<K; jj++)
			{
#pragma omp simd
        for(std::size_t j=0; j<M; j++)
					tmpJ[i][j] += a[jj][j] * b[jj][i];
			}
		}
	

		// Check if the interpolated solution is needed.
		if( C )
		{
			// Cast array C from 1D to 2D.
			as3double (*c)[M] = (as3double (*)[M]) &C[l*M2];

			for(std::size_t j=0; j<M; j++)
			{
        for(std::size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
          for(std::size_t i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to C.
			for(std::size_t j=0; j<M; j++)
				for(std::size_t i=0; i<M; i++)
					c[j][i] = tmpI[j][i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Cast array CDerR from 1D to 2D.
			as3double (*cDerR)[M] = (as3double (*)[M]) &CDerR[l*M2];

			for(std::size_t j=0; j<M; j++)
			{
				for(std::size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
          for(std::size_t i=0; i<M; i++)
						tmpI[j][i] += aDer[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerR.
			for(std::size_t j=0; j<M; j++)
				for(std::size_t i=0; i<M; i++)
					cDerR[j][i] = tmpI[j][i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Cast array CDerS from 1D to 2D. 
			as3double (*cDerS)[M] = (as3double (*)[M]) &CDerS[l*M2];

			for(std::size_t i=0; i<K; i++)
			{
        for(std::size_t s=0; s<M; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<K; jj++)
				{
#pragma omp simd
          for(std::size_t j=0; j<M; j++)
						tmpJ[i][j] += aDer[jj][j] * b[jj][i];
				}
			}

			for(std::size_t j=0; j<M; j++)
			{
        for(std::size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
          for(std::size_t i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerS.
			for(std::size_t j=0; j<M; j++)
				for(std::size_t i=0; i<M; i++)
					cDerS[j][i] = tmpI[j][i];
		}
	}
}



