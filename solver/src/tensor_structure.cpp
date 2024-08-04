#include "tensor_structure.hpp"


//-----------------------------------------------------------------------------------
// ITensorProduct member functions.
//-----------------------------------------------------------------------------------


void ITensorProduct::CustomVolume
(
 size_t           K,
 size_t           N,
 size_t           M,
 const as3double *A,
 const as3double *ADer,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that implements a generic (run-time) tensor-product volume implementation.
	*/
{
	// Total number of (volume) input/output points in 2D.
	const size_t M2 = M*M; 
	const size_t K2 = K*K;

	
	// Loop over each N entry.
	for(size_t l=0; l<N; l++)
	{
		// Create temporary storage.
		as3double tmpJ[K][M];
		as3double tmpI[M][M];

		// Cast array B from 1D to 2D.
		const as3double *b = &B[l*K2];

		for(size_t i=0; i<K; i++)
		{
			for(size_t s=0; s<M; s++) tmpJ[i][s] = 0.0;
			for(size_t jj=0; jj<K; jj++)
			{
        for(size_t j=0; j<M; j++)
					tmpJ[i][j] += A[jj*M+j] * b[jj*K+i];
			}
		}
	

		// Check if the interpolated solution is needed.
		if( C )
		{
			// Cast array C from 1D to 2D.
			as3double *c = &C[l*M2];

			for(size_t j=0; j<M; j++)
			{
        for(size_t s=0; s<M; s++) tmpI[j][s] = 0.0;
				for(size_t ii=0; ii<K; ii++)
				{
          for(size_t i=0; i<M; i++)
						tmpI[j][i] += A[ii*M+i] * tmpJ[ii][j]; 
				}
			}

			// Copy values to C.
			for(size_t j=0; j<M; j++)
				for(size_t i=0; i<M; i++)
					c[j*M+i] = tmpI[j][i]; 
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Cast array CDerR from 1D to 2D.
			as3double *cDerR = &CDerR[l*M2]; 

			for(size_t j=0; j<M; j++)
			{
				for(size_t s=0; s<M; s++) tmpI[j][s] = 0.0;
				for(size_t ii=0; ii<K; ii++)
				{
          for(size_t i=0; i<M; i++)
						tmpI[j][i] += ADer[ii*M+i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerR.
			for(size_t j=0; j<M; j++)
				for(size_t i=0; i<M; i++)
					cDerR[j*M+i] = tmpI[j][i];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Cast array CDerS from 1D to 2D. 
			as3double *cDerS = &CDerS[l*M2];

			for(size_t i=0; i<K; i++)
			{
        for(size_t s=0; s<M; s++) tmpJ[i][s] = 0.0;
				for(size_t jj=0; jj<K; jj++)
				{
          for(size_t j=0; j<M; j++)
						tmpJ[i][j] += ADer[jj*M+j] * b[jj*K+i];
				}
			}

			for(size_t j=0; j<M; j++)
			{
        for(size_t s=0; s<M; s++) tmpI[j][s] = 0.0;
				for(size_t ii=0; ii<K; ii++)
				{
          for(size_t i=0; i<M; i++)
						tmpI[j][i] += A[ii*M+i] * tmpJ[ii][j];
				}
			}

			// Copy values to CDerS.
			for(size_t j=0; j<M; j++)
				for(size_t i=0; i<M; i++)
					cDerS[j*M+i] = tmpI[j][i];
		}
	}
}

//-----------------------------------------------------------------------------------

void ITensorProduct::CustomSurfaceIMIN
(
 size_t           K,
 size_t           N,
 size_t           M,
 const as3double *A,
 const as3double *ADer,
 const as3double *ADerFace,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that implements a generic (run-time) tensor-product IMIN surface implementation.
	*/
{
	// Total number of (volume) input points in 2D.
	const size_t K2 = K*K;


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

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
				for(size_t jj=0; jj<K; jj++)
					tmpI[k] += A[jj*M+k] * tmpJ[jj];
			}

			// Copy values to C.
      for(size_t j=0; j<M; j++) c[j] = tmpI[j];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpI[k] += ADer[jj*M+k] * tmpJ[jj];
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
			const as3double *b = &B[l*K2];

			for(size_t s=0; s<K; s++) tmpJ[s] = 0.0;
			for(size_t j=0; j<K; j++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpJ[j] += ADerFace[ii] * b[j*K+ii];
			}

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpI[k] += A[jj*M+k] * tmpJ[jj];
			}

			// Copy values to CDerR.
			for(size_t j=0; j<M; j++) cDerR[j] = tmpI[j];
		}
	}
}

//-----------------------------------------------------------------------------------

void ITensorProduct::CustomSurfaceIMAX
(
 size_t           K,
 size_t           N,
 size_t           M,
 const as3double *A,
 const as3double *ADer,
 const as3double *ADerFace,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that implements a generic (run-time) tensor-product IMAX surface implementation.
	*/
{
	// Total number of (volume) input points in 2D.
	const size_t K2  = K*K;
	const size_t Km1 = K-1;

	// Loop over each N entry.
	for(size_t l=0; l<N; l++)
	{
		// Create temporary storage.
		as3double tmpI[M];
		as3double tmpJ[K];

		// Store the IMIN boundary surface data in tmpJ.
    for(size_t s=0; s<K; s++) tmpJ[s] = B[l*K2+K*s+Km1];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
				for(size_t jj=0; jj<K; jj++)
					tmpI[k] += A[jj*M+k] * tmpJ[jj];
			}

			// Copy values to C.
      for(size_t j=0; j<M; j++) c[j] = tmpI[j];
		}


		// Check if the derivative in the s-direction is needed.
		if( CDerS )
		{
			// Get a pointer to the current variable address in CDerS.
			as3double *cDerS = &CDerS[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpI[k] += ADer[jj*M+k] * tmpJ[jj];
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
			const as3double *b = &B[l*K2];

			for(size_t s=0; s<K; s++) tmpJ[s] = 0.0;
			for(size_t j=0; j<K; j++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpJ[j] += ADerFace[ii] * b[j*K+ii];
			}

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpI[k] += A[jj*M+k] * tmpJ[jj];
			}

			// Copy values to CDerR.
			for(size_t j=0; j<M; j++) cDerR[j] = tmpI[j];
		}
	}
}

//-----------------------------------------------------------------------------------

void ITensorProduct::CustomSurfaceJMIN
(
 size_t           K,
 size_t           N,
 size_t           M,
 const as3double *A,
 const as3double *ADer,
 const as3double *ADerFace,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that implements a generic (run-time) tensor-product JMIN surface implementation.
	*/
{
	// Total number of (volume) solution points in 2D.
	const size_t K2 = K*K;


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

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += A[ii*M+k] * tmpJ[ii];
			}

			// Copy values to C.
			for(size_t i=0; i<M; i++) c[i] = tmpI[i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += ADer[ii*M+k] * tmpJ[ii];
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
			const as3double *b = &B[l*K2];

			for(size_t s=0; s<K; s++) tmpJ[s] = 0.0;
			for(size_t i=0; i<K; i++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpJ[i] += ADerFace[jj] * b[jj*K+i];
			}

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += A[ii*M+k] * tmpJ[ii];
			}

			// Copy values to CDerS.
			for(size_t i=0; i<M; i++) cDerS[i] = tmpI[i];
		}
	}
}

//-----------------------------------------------------------------------------------

void ITensorProduct::CustomSurfaceJMAX
(
 size_t           K,
 size_t           N,
 size_t           M,
 const as3double *A,
 const as3double *ADer,
 const as3double *ADerFace,
 const as3double *B,
 as3double       *C,
 as3double       *CDerR,
 as3double       *CDerS
)
 /*
	* Function that implements a generic (run-time) tensor-product JMAX surface implementation.
	*/
{
	// Total number of (volume) solution points in 2D.
	const size_t K2   = K*K;
	const size_t K2mK = K2-K;


	// Loop over each N entry.
	for(size_t l=0; l<N; l++)
	{
		// Create temporary storage.
		as3double tmpI[M];
		as3double tmpJ[K];

		// Store boundary surface info in tmpJ.
    for(size_t s=0; s<K; s++) tmpJ[s] = B[l*K2+s+K2mK];


		// Check if the interpolated solution is needed.
		if( C )
		{
			// Get a pointer to the current variable address in C.
			as3double *c = &C[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += A[ii*M+k] * tmpJ[ii];
			}

			// Copy values to C.
			for(size_t i=0; i<M; i++) c[i] = tmpI[i];
		}


		// Check if the derivative in the r-direction is needed.
		if( CDerR )
		{
			// Get a pointer to the current variable address in CDerR.
			as3double *cDerR = &CDerR[l*M];

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += ADer[ii*M+k] * tmpJ[ii];
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
			const as3double *b = &B[l*K2];

			for(size_t s=0; s<K; s++) tmpJ[s] = 0.0;
			for(size_t i=0; i<K; i++)
			{
        for(size_t jj=0; jj<K; jj++)
					tmpJ[i] += ADerFace[jj] * b[jj*K+i];
			}

			for(size_t s=0; s<M; s++) tmpI[s] = 0.0;
			for(size_t k=0; k<M; k++)
			{
        for(size_t ii=0; ii<K; ii++)
					tmpI[k] += A[ii*M+k] * tmpJ[ii];
			}

			// Copy values to CDerS.
			for(size_t i=0; i<M; i++) cDerS[i] = tmpI[i];
		}
	}
}


