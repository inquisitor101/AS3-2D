/*!
 * @file
 * @brief Collocated compile-time tensor kernels for volume.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the volume templated functions in CCollocatedTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the volume integration
 *   points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B solution values [N][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [N][K*K]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [N][K*K]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [N][K*K]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::CompileTimeVolume
(
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that evaluates collocated values and reference derivatives in the volume.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *ADer = mDerLagrangeInt1D.data();

	// Cast the differentiation operator from 1D to 2D.
	const as3double (*aDer)[K] = (const as3double (*)[K]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr std::size_t K2 = K*K;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Cast array C from 1D to 2D.
			as3double (*c)[K] = (as3double (*)[K]) &C[l*K2];

			// Collocated interpolation is the identity.
			for(std::size_t j=0; j<K; j++)
			{
#pragma omp simd
				for(std::size_t i=0; i<K; i++)
					c[j][i] = b[j][i];
			}
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Cast array CDerR from 1D to 2D.
			as3double (*cDerR)[K] = (as3double (*)[K]) &CDerR[l*K2];

			// Tensor product in the r-direction; the s-interpolation is the identity.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t s=0; s<K; s++) cDerR[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
					for(std::size_t i=0; i<K; i++)
						cDerR[j][i] += aDer[i][ii] * b[j][ii];
				}
			}
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Cast array CDerS from 1D to 2D.
			as3double (*cDerS)[K] = (as3double (*)[K]) &CDerS[l*K2];

			// Tensor product in the s-direction; the r-interpolation is the identity.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t s=0; s<K; s++) cDerS[j][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<K; jj++)
				{
#pragma omp simd
					for(std::size_t i=0; i<K; i++)
						cDerS[j][i] += aDer[j][jj] * b[jj][i];
				}
			}
		}
	}
}


//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the volume integration
 *   points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] n number of variables processed with the fixed K operators.
 * @param[in] B solution values [n][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [n][K*K]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [n][K*K]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [n][K*K]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::PartialCompileTimeVolume
(
 std::size_t      n,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that evaluates the volume operators for the requested number of variables.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *ADer = mDerLagrangeInt1D.data();

	// Cast the differentiation operator from 1D to 2D.
	const as3double (*aDer)[K] = (const as3double (*)[K]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr std::size_t K2 = K*K;


	// Loop over each n entry.
	for(std::size_t l=0; l<n; l++)
	{
		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Cast array C from 1D to 2D.
			as3double (*c)[K] = (as3double (*)[K]) &C[l*K2];

			// Collocated interpolation is the identity.
			for(std::size_t j=0; j<K; j++)
			{
#pragma omp simd
				for(std::size_t i=0; i<K; i++)
					c[j][i] = b[j][i];
			}
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Cast array CDerR from 1D to 2D.
			as3double (*cDerR)[K] = (as3double (*)[K]) &CDerR[l*K2];

			// Tensor product in the r-direction; the s-interpolation is the identity.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t s=0; s<K; s++) cDerR[j][s] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
					for(std::size_t i=0; i<K; i++)
						cDerR[j][i] += aDer[i][ii] * b[j][ii];
				}
			}
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Cast array CDerS from 1D to 2D.
			as3double (*cDerS)[K] = (as3double (*)[K]) &CDerS[l*K2];

			// Tensor product in the s-direction; the r-interpolation is the identity.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t s=0; s<K; s++) cDerS[j][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<K; jj++)
				{
#pragma omp simd
					for(std::size_t i=0; i<K; i++)
						cDerS[j][i] += aDer[j][jj] * b[jj][i];
				}
			}
		}
	}
}
