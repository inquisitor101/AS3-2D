/*!
 * @file
 * @brief Collocated compile-time tensor kernels for residual surface imax.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the residual surface imax templated functions in CCollocatedTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that lifts weighted inputs from the IMAX surface integration points to the solution DOFs.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] B optional weighted value contribution [N][K]; nullptr omits this contribution.
 * @param[in] BDerR optional weighted r-derivative contribution [N][K]; nullptr omits this contribution.
 * @param[in] BDerS optional weighted s-derivative contribution [N][K]; nullptr omits this contribution.
 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
 *
 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::CompileTimeResidualSurfaceIMAX
(
 const as3double *B,
 const as3double *BDerR,
 const as3double *BDerS,
 as3double       *C
)
 /*
	* Function that accumulates weighted IMAX face contributions into the residual.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *ADer = mDerLagrangeInt1D.data();

	// Cast the differentiation operator from 1D to 2D.
	const as3double (*aDer)[K] = (const as3double (*)[K]) ADer;

	// Total number of (volume) solution points in 2D and index of the IMAX node.
	constexpr std::size_t K2 = K*K;
	constexpr std::size_t IK = K-1;


	// Loop over each N entry.
	for(std::size_t l=0; l<N; l++)
	{
		// Cast array C from 1D to 2D.
		as3double (*c)[K] = (as3double (*)[K]) &C[l*K2];


		// Check if the residual includes a source term.
		if( B )
		{
			// Get a pointer to the current variable address in B.
			const as3double *b = &B[l*K];

#pragma omp simd
			for(std::size_t k=0; k<K; k++)
				c[k][IK] += b[k];
		}


		// Check if the residual includes an r-differential term.
		if( BDerR )
		{
			// Get a pointer to the current variable address in BDerR.
			const as3double *bDerR = &BDerR[l*K];

			// Distribute the endpoint derivative over the volume solution nodes.
			for(std::size_t j=0; j<K; j++)
			{
#pragma omp simd
				for(std::size_t i=0; i<K; i++)
					c[j][i] += aDer[IK][i] * bDerR[j];
			}
		}


		// Check if the residual includes an s-differential term.
		if( BDerS )
		{
			// Get a pointer to the current variable address in BDerS.
			const as3double *bDerS = &BDerS[l*K];

			// Apply the transposed derivative along the face.
#pragma omp simd
			for(std::size_t k=0; k<K; k++)
			{
				as3double value = static_cast<as3double>(0.0);
				for(std::size_t ii=0; ii<K; ii++)
					value += aDer[ii][k] * bDerS[ii];
				c[k][IK] += value;
			}
		}
	}
}
