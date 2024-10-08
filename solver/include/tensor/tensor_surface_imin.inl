

//-----------------------------------------------------------------------------------
// Implementation of the IMIN surface templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


template<size_t K, size_t M>
void CTensorProduct<K, M>::SurfaceIMIN
(
 const size_t     N,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that computes the tensor product on the surface integration points of IMIN.
	*/
{
	// Unpack the information needed, based on the standard element.
	const as3double *A        = mLagrangeInt1DTrans.data();
	const as3double *ADer     = mDerLagrangeInt1DTrans.data();
	const as3double *aDerFace = mDerLagrangeMinFace1D.data();

	// Cast the arrays from 1D to 2D.
	const as3double (*a)[M]    = (const as3double (*)[M]) A;
	const as3double (*aDer)[M] = (const as3double (*)[M]) ADer;

	// Total number of (volume) solution points in 2D.
	constexpr size_t K2 = K*K;


	// Loop over each N entry.
	for(size_t l=0; l<N; l++)
	{
		// Create temporary storage.
		as3double tmpI[M];
		as3double tmpJ[K];

		// Store the IMIN boundary surface data in tmpJ.
    for(size_t s=0; s<K; s++) tmpJ[s] = B[l*K2+s*K];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
			for(size_t k=0; k<M; k++)
			{
				for(size_t jj=0; jj<K; jj++)
					tmpI[k] += a[jj][k] * tmpJ[jj];
			}

			// Copy values to C.
      for(size_t j=0; j<M; j++) c[j] = tmpI[j];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
			for(size_t k=0; k<M; k++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpI[k] += aDer[jj][k] * tmpJ[jj];
			}

			// Copy values to CDerS.
      for(size_t j=0; j<M; j++) cDerS[j] = tmpI[j];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*M];

			// Cast array B from 1D to 2D.
			const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];

			for(size_t s=0; s<K; s++) tmpJ[s] = static_cast<as3double>(0.0);
			for(size_t j=0; j<K; j++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpJ[j] += aDerFace[ii] * b[j][ii];
			}

			for(size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
			for(size_t k=0; k<M; k++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpI[k] += a[jj][k] * tmpJ[jj];
			}

			// Copy values to CDerR.
			for(size_t j=0; j<M; j++) cDerR[j] = tmpI[j];
		}
	}
}
