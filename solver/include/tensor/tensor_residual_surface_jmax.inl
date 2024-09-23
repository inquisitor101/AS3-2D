

//-----------------------------------------------------------------------------------
// Implementation of the JMAX surface residual templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


template<size_t K, size_t M>
void CTensorProduct<K, M>::ResidualSurfaceJMAX
(
 const size_t     N,
 const as3double *B,
 const as3double *BDerR,
 const as3double *BDerS,
 as3double       *C
)
 /*
	* Function that computes the residual via a tensor product on the JMAX surface integration points.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A        = mLagrangeInt1DTrans.data();
	const as3double *ADer     = mDerLagrangeInt1DTrans.data();
	const as3double *aDerFace = mDerLagrangeMaxFace1D.data();

	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) integration/solution points in 2D.
	constexpr size_t K2   = K*K;
	// Index of the JMAX node.
	constexpr size_t JK = K-1;

	// Ensure the tensor dimension is correct.
#if DEBUG	
	if( (K != mStandardElementContainer->GetnSol1D()) 
			|| 
			(M != mStandardElementContainer->GetnInt1D()) ) 
		ERROR("Mismatch in tensor dimensions.");
#endif


	// Loop over each N entry.
	for(size_t l=0; l<N; l++)
	{
    // Cast array C from 1D to 2D.
    as3double (*c)[K]  = (as3double (*)[K]) &C[l*K2];


		// Check if the residual includes a source term.
		if( B )
		{
			// Get a pointer to the current variable address in B.
			const as3double *b = &B[l*M];
    	
			for(size_t k=0; k<K; k++)
			{
    	  for(size_t ii=0; ii<M; ii++)
    	    c[JK][k] += a[k][ii] * b[ii];
    	}
		}


		// Check if the residual includes an r-differential term.
		if( BDerR )
		{
			// Get a pointer to the current variable address in BDerR.
			const as3double *bDerR = &BDerR[l*M];
  
    	for(size_t k=0; k<K; k++)
			{
    	  for(size_t ii=0; ii<M; ii++)
    	    c[JK][k] += aDer[k][ii] * bDerR[ii];
    	}
		}

		// Check if the residual includes an s-differential term.
		if( BDerS )
		{
			// Get a pointer to the current variable address in BDerS.
			const as3double *bDerS = &BDerS[l*M];

			// Create and initialize a temporary storage to zero.
			as3double tmpI[K] = {};

			// Interpolate the solution to the r-nodes.
			for(size_t s=0; s<K; s++)
			{
				for(size_t jj=0; jj<M; jj++)
					tmpI[s] += a[s][jj] * bDerS[jj]; 
			}

			// Update the residual, by including the derivative in the s-direction.
			for(size_t i=0; i<K; i++)
				for(size_t j=0; j<K; j++)
					c[i][j] += aDerFace[i] * tmpI[j];

		}
	}
}
