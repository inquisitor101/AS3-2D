

//-----------------------------------------------------------------------------------
// Implementation of the JMIN surface templated function in CTensorProduct.
//-----------------------------------------------------------------------------------


template<size_t K, size_t M>
void CTensorProduct<K, M>::SurfaceJMIN
(
 const size_t     N,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that computes the tensor product on the surface integration points of JMIN.
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
		as3double tmpI[M];
		as3double tmpJ[K];

		// Store boundary surface info in tmpJ.
    for(size_t s=0; s<K; s++) tmpJ[s] = B[l*K2+s];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to C.
			for(size_t i=0; i<M; i++) c[i] = tmpI[i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += aDer[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerR.
			for(size_t i=0; i<M; i++) cDerR[i] = tmpI[i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*M];

			// Cast array B from 1D to 2D.
			const as3double (*b)[K] = (const as3double (*)[K]) &B[l*K2];
			
			for(size_t s=0; s<K; s++) tmpJ[s] = static_cast<as3double>(0.0);
			for(size_t i=0; i<K; i++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpJ[i] += aDerFace[jj] * b[jj][i];
			}

			for(size_t s=0; s<M; s++) tmpI[s] = static_cast<as3double>(0.0);
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += a[ii][k] * tmpJ[ii];
			}

			// Copy values to CDerS.
			for(size_t i=0; i<M; i++) cDerS[i] = tmpI[i];
		}
	}
}
