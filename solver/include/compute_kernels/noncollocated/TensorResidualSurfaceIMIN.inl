/*!
 * @file
 * @brief Compile-time tensor kernels for residual surface imin.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the IMIN surface residual templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that lifts weighted inputs from the IMIN surface integration points to the solution DOFs.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this contribution.
 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this contribution.
 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::CompileTimeResidualSurfaceIMIN
(
 const as3double *B,
 const as3double *BDerR,
 const as3double *BDerS,
 as3double       *C
)
 /*
	* Function that computes the residual via a tensor product on the IMIN surface integration points.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A        = mLagrangeInt1DTrans.data();
	const as3double *ADer     = mDerLagrangeInt1DTrans.data();
	const as3double *aDerFace = mDerLagrangeMinFace1D.data();

	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr std::size_t K2 = K*K;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
    // Cast array C from 1D to 2D.
    as3double (*c)[K]  = (as3double (*)[K]) &C[l*K2];


		// Check if the residual includes a source term.
		if( B )
		{
			// Get a pointer to the current variable address in B.
			const as3double *b = &B[l*M];
    	
			for(std::size_t k=0; k<K; k++)
			{
    	  for(std::size_t ii=0; ii<M; ii++)
    	    c[k][0] += a[k][ii] * b[ii];
    	}
		}


		// Check if the residual includes an r-differential term.
		if( BDerR )
		{
			// Get a pointer to the current variable address in BDerR.
			const as3double *bDerR = &BDerR[l*M];

			// Create and initialize a temporary storage to zero.
			as3double tmpJ[K] = {};

			// Interpolate the solution to the s-nodes.
			for(std::size_t s=0; s<K; s++)
			{
				for(std::size_t jj=0; jj<M; jj++)
					tmpJ[s] += a[s][jj] * bDerR[jj]; 
			}

			// Update the residual, by including the derivative in the r-direction.
			for(std::size_t i=0; i<K; i++)
			{
#pragma omp simd
				for(std::size_t j=0; j<K; j++)
					c[i][j] += aDerFace[j] * tmpJ[i];
			}
		}


		// Check if the residual includes an s-differential term.
		if( BDerS )
		{
			// Get a pointer to the current variable address in BDerS.
			const as3double *bDerS = &BDerS[l*M];

    	for(std::size_t k=0; k<K; k++)
			{
    	  for(std::size_t ii=0; ii<M; ii++)
    	    c[k][0] += aDer[k][ii] * bDerS[ii];
    	}
		}
	}
}
