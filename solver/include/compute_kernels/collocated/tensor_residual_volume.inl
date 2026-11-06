/*!
 * @file
 * @brief Collocated compile-time tensor kernels for residual volume.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the residual volume templated functions in CCollocatedTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that lifts weighted inputs from the volume integration points to the solution DOFs.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B optional weighted value contribution [N][K*K]; nullptr omits this contribution.
 * @param[in] BDerR optional weighted r-derivative contribution [N][K*K]; nullptr omits this contribution.
 * @param[in] BDerS optional weighted s-derivative contribution [N][K*K]; nullptr omits this contribution.
 * @param[out] C solution residual [N][K*K], overwritten by the volume contribution.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::CompileTimeResidualVolume
(
 const as3double *B,
 const as3double *BDerR,
 const as3double *BDerS,
 as3double       *C
)
 /*
	* Function that computes the collocated volume residual and resets the output.
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
		// Cast array C from 1D to 2D.
		as3double (*c)[K] = (as3double (*)[K]) &C[l*K2];

		// Reset the residual before adding the volume contributions.
		for(std::size_t j=0; j<K; j++)
		{
#pragma omp simd
			for(std::size_t i=0; i<K; i++)
				c[j][i] = static_cast<as3double>(0.0);
		}


		// Check if the residual includes a source term.
		if( B )
		{
			// Cast array B from 1D to 2D.
			const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];

			// Collocated value lifting is the identity.
			for(std::size_t j=0; j<K; j++)
			{
#pragma omp simd
				for(std::size_t i=0; i<K; i++)
					c[j][i] = b[j][i];
			}
		}


		// Check if the residual includes an r-differential term.
		if( BDerR )
		{
			// Cast array BDerR from 1D to 2D.
			const as3double (*bDerR)[K] = (const as3double (*)[K]) &BDerR[l*K2];

			// Transposed tensor product in the r-direction.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t ii=0; ii<K; ii++)
				{
#pragma omp simd
					for(std::size_t i=0; i<K; i++)
						c[j][i] += aDer[ii][i] * bDerR[j][ii];
				}
			}
		}


		// Check if the residual includes an s-differential term.
		if( BDerS )
		{
			// Cast array BDerS from 1D to 2D.
			const as3double (*bDerS)[K] = (const as3double (*)[K]) &BDerS[l*K2];

			// Transposed tensor product in the s-direction.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t jj=0; jj<K; jj++)
				{
#pragma omp simd
					for(std::size_t i=0; i<K; i++)
						c[j][i] += aDer[jj][j] * bDerS[jj][i];
				}
			}
		}
	}
}
