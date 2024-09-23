

//-----------------------------------------------------------------------------------
// Implementation of the volume templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


template<size_t K, size_t M>
void CTensorProduct<K, M>::Volume
(
 const size_t     N,
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
	constexpr size_t M2 = M*M; 
	constexpr size_t K2 = K*K;

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
		// Create temporary storage.
		as3double tmpJ[K][M];
		as3double tmpI[M][M];

		// Cast array B from 1D to 2D.
		const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];

		for(size_t i=0; i<K; i++)
		{
			for(size_t s=0; s<M; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
			for(size_t jj=0; jj<K; jj++)
			{
        for(size_t j=0; j<M; j++)
					tmpJ[i][j] += a[jj][j] * b[jj][i];
			}
		}
	

		// Check if the interpolated solution is needed.
		if( C )
		{
			// Cast array C from 1D to 2D.
			as3double (*c)[M] = (as3double (*)[M]) &C[l*M2];

			for(size_t j=0; j<M; j++)
			{
        for(size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(size_t ii=0; ii<K; ii++)
				{
          for(size_t i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to C.
			for(size_t j=0; j<M; j++)
				for(size_t i=0; i<M; i++)
					c[j][i] = tmpI[j][i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Cast array CDerR from 1D to 2D.
			as3double (*cDerR)[M] = (as3double (*)[M]) &CDerR[l*M2];

			for(size_t j=0; j<M; j++)
			{
				for(size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(size_t ii=0; ii<K; ii++)
				{
          for(size_t i=0; i<M; i++)
						tmpI[j][i] += aDer[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerR.
			for(size_t j=0; j<M; j++)
				for(size_t i=0; i<M; i++)
					cDerR[j][i] = tmpI[j][i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Cast array CDerS from 1D to 2D. 
			as3double (*cDerS)[M] = (as3double (*)[M]) &CDerS[l*M2];

			for(size_t i=0; i<K; i++)
			{
        for(size_t s=0; s<M; s++) tmpJ[i][s] = static_cast<as3double>(0.0);
				for(size_t jj=0; jj<K; jj++)
				{
          for(size_t j=0; j<M; j++)
						tmpJ[i][j] += aDer[jj][j] * b[jj][i];
				}
			}

			for(size_t j=0; j<M; j++)
			{
        for(size_t s=0; s<M; s++) tmpI[j][s] = static_cast<as3double>(0.0);
				for(size_t ii=0; ii<K; ii++)
				{
          for(size_t i=0; i<M; i++)
						tmpI[j][i] += a[ii][i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerS.
			for(size_t j=0; j<M; j++)
				for(size_t i=0; i<M; i++)
					cDerS[j][i] = tmpI[j][i];
		}
	}
}



