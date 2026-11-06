/*!
 * @file
 * @brief Collocated compile-time tensor kernels for surface jmin.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the surface jmin templated functions in CCollocatedTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the JMIN surface
 *   integration points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B solution values [N][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [N][K]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [N][K]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [N][K]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::CompileTimeSurfaceJMIN
(
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that evaluates values and reference derivatives on the JMIN face.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *ADer = mDerLagrangeInt1D.data();

	// Cast the differentiation operator from 1D to 2D.
	const as3double (*aDer)[K] = (const as3double (*)[K]) ADer;

	// Total number of (volume) solution points in 2D and index of the JMIN node.
	constexpr std::size_t K2 = K*K;
	constexpr std::size_t JK = 0;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*K];

			// Collocated face interpolation is a direct boundary gather.
#pragma omp simd
			for(std::size_t k=0; k<K; k++)
				c[k] = b[JK][k];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*K];

#pragma omp simd
			for(std::size_t k=0; k<K; k++)
			{
				cDerR[k] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
					cDerR[k] += aDer[k][ii] * b[JK][ii];
			}
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*K];

#pragma omp simd
			for(std::size_t k=0; k<K; k++)
			{
				cDerS[k] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<K; jj++)
					cDerS[k] += aDer[JK][jj] * b[jj][k];
			}
		}
	}
}


//-----------------------------------------------------------------------------------


/*!
 * @brief Function that interpolates the solution and optional parametric derivatives to the JMIN surface
 *   integration points.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] n number of variables processed with the fixed K operators.
 * @param[in] B solution values [n][K*K], with i varying fastest.
 * @param[out] C optional interpolated values [n][K]; nullptr skips this output.
 * @param[out] CDerR optional derivative in the r-direction [n][K]; nullptr skips this output.
 * @param[out] CDerS optional derivative in the s-direction [n][K]; nullptr skips this output.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::PartialCompileTimeSurfaceJMIN
(
 std::size_t      n,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that evaluates the JMIN face operators for the requested number of variables.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *ADer = mDerLagrangeInt1D.data();

	// Cast the differentiation operator from 1D to 2D.
	const as3double (*aDer)[K] = (const as3double (*)[K]) ADer;

	// Total number of (volume) solution points in 2D and index of the JMIN node.
	constexpr std::size_t K2 = K*K;
	constexpr std::size_t JK = 0;


	// Loop over each n entry.
	for(std::size_t l=0; l<n; l++)
	{
		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*K];

			// Collocated face interpolation is a direct boundary gather.
#pragma omp simd
			for(std::size_t k=0; k<K; k++)
				c[k] = b[JK][k];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*K];

#pragma omp simd
			for(std::size_t k=0; k<K; k++)
			{
				cDerR[k] = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
					cDerR[k] += aDer[k][ii] * b[JK][ii];
			}
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*K];

#pragma omp simd
			for(std::size_t k=0; k<K; k++)
			{
				cDerS[k] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<K; jj++)
					cDerS[k] += aDer[JK][jj] * b[jj][k];
			}
		}
	}
}
