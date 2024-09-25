

//-----------------------------------------------------------------------------------
// Implementation of the volume residual templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


template<size_t K, size_t M>
void CTensorProduct<K, M>::ResidualVolume
(
 const size_t     N,
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
	constexpr size_t M2 = M*M; 
	constexpr size_t K2 = K*K;


	// Loop over each N entry.
	for(size_t l=0; l<N; l++)
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
			for(size_t i=0; i<M; i++)
			{
  			for(size_t s=0; s<K; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(size_t jj=0; jj<M; jj++)
				{
    	    for(size_t j=0; j<K; j++)
						tmpJ[i][j] += a[jj][j] * b[jj][i];
				}
			}
	
			// Tensor product in the r-direction.
			for(size_t j=0; j<K; j++)
			{
				for(size_t ii=0; ii<M; ii++)
				{
    	    for(size_t i=0; i<K; i++)
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
			for(size_t i=0; i<M; i++)
			{
  			for(size_t s=0; s<K; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(size_t jj=0; jj<M; jj++)
				{
    	    for(size_t j=0; j<K; j++)
						tmpJ[i][j] += a[jj][j] * bDerR[jj][i];
				}
			}

			// Tensor product in the r-direction.
			for(size_t j=0; j<K; j++)
			{
				for(size_t ii=0; ii<M; ii++)
				{
    	    for(size_t i=0; i<K; i++)
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
			for(size_t i=0; i<M; i++)
			{
  			for(size_t s=0; s<K; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(size_t jj=0; jj<M; jj++)
				{
    	    for(size_t j=0; j<K; j++)
						tmpJ[i][j] += aDer[jj][j] * bDerS[jj][i];
				}
			}

			// Tensor product in the r-direction.
			for(size_t j=0; j<K; j++)
			{
				for(size_t ii=0; ii<M; ii++)
				{
    	    for(size_t i=0; i<K; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}
		}

		// Copy values to C. Notice, this overrides the existing data in C.
		for(unsigned short j=0; j<K; j++)
			for(unsigned short i=0; i<K; i++)
				c[j][i] = tmpI[j][i];
	}
}
