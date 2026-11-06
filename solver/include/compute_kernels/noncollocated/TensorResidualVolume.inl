/*!
 * @file
 * @brief Compile-time tensor kernels for residual volume.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the volume residual templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that lifts weighted inputs from the volume integration points to the solution DOFs.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B optional weighted value contribution [N][M*M]; nullptr omits this contribution.
 * @param[in] BDerR optional weighted r-derivative contribution [N][M*M]; nullptr omits this contribution.
 * @param[in] BDerS optional weighted s-derivative contribution [N][M*M]; nullptr omits this contribution.
 * @param[out] C solution residual [N][K*K], overwritten by the volume contribution.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::CompileTimeResidualVolume
(
 const as3double *B,
 const as3double *BDerR,
 const as3double *BDerS,
 as3double       *C
)
 /*
	* Function that computes the residual via a tensor product on the volume integration points.
	* NOTE, this also resets the residual.
	*/
{
 	// Unpack the information needed, based on the standard element.
	const as3double *A    = mLagrangeInt1D.data();
	const as3double *ADer = mDerLagrangeInt1D.data();

	// Cast the arrays from 1D to 2D, for convenience.
	const as3double (*a)[K]    = (const as3double (*)[K]) A;
	const as3double (*aDer)[K] = (const as3double (*)[K]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr std::size_t M2 = M*M; 
	constexpr std::size_t K2 = K*K;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
		// Create temporary storage and ensure tmpI is initialized to zero.
		as3double tmpI[K][K] = {}; 
		as3double tmpJ[M][K];

    // Cast array C from 1D to 2D.
    as3double (*c)[K] = (as3double (*)[K]) &C[l*K2];

		
		// Check if a source term is required.
		if( B )
		{
			// Cast array B from 1D to 2D. 
			const as3double (*b)[M] = (const as3double (*)[M]) &B[l*M2];

			// Tensor product in the s-direction.
			for(std::size_t i=0; i<M; i++)
			{
  			for(std::size_t s=0; s<K; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<M; jj++)
				{
#pragma omp simd
    	    for(std::size_t j=0; j<K; j++)
						tmpJ[i][j] += a[jj][j] * b[jj][i];
				}
			}
	
			// Tensor product in the r-direction.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t ii=0; ii<M; ii++)
				{
#pragma omp simd
    	    for(std::size_t i=0; i<K; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}
		}


		// Check if an r-differential term is required.
		if( BDerR )
		{
			// Cast array BDerR from 1D to 2D. 
			const as3double (*bDerR)[M] = (const as3double (*)[M]) &BDerR[l*M2];

			// Tensor product in the s-direction.
			for(std::size_t i=0; i<M; i++)
			{
  			for(std::size_t s=0; s<K; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<M; jj++)
				{
#pragma omp simd
    	    for(std::size_t j=0; j<K; j++)
						tmpJ[i][j] += a[jj][j] * bDerR[jj][i];
				}
			}

			// Tensor product in the r-direction.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t ii=0; ii<M; ii++)
				{
#pragma omp simd
    	    for(std::size_t i=0; i<K; i++)
						tmpI[j][i] += aDer[ii][i] * tmpJ[ii][j];
				}
			}
		}


		// Check if an s-differential term is required.
		if( BDerS )
		{
			// Cast array BDerS from 1D to 2D. 
			const as3double (*bDerS)[M] = (const as3double (*)[M]) &BDerS[l*M2];

			// Tensor product in the s-direction.
			for(std::size_t i=0; i<M; i++)
			{
  			for(std::size_t s=0; s<K; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(std::size_t jj=0; jj<M; jj++)
				{
#pragma omp simd
    	    for(std::size_t j=0; j<K; j++)
						tmpJ[i][j] += aDer[jj][j] * bDerS[jj][i];
				}
			}

			// Tensor product in the r-direction.
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t ii=0; ii<M; ii++)
				{
#pragma omp simd
    	    for(std::size_t i=0; i<K; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}
		}

		// Copy values to C. Notice, this overrides the existing data in C.
		for(std::size_t j=0; j<K; j++)
			for(std::size_t i=0; i<K; i++)
				c[j][i] = tmpI[j][i];
	}
}
