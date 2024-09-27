#include "riemann_solver_structure.hpp"



//-----------------------------------------------------------------------------------
// IRiemannSolver member functions.
//-----------------------------------------------------------------------------------


IRiemannSolver::IRiemannSolver
(
 CConfig *config_container
)
 /*
	* Constructor for the interface Riemann solver class.
	*/
{

}

//-----------------------------------------------------------------------------------

IRiemannSolver::~IRiemannSolver
(
 void
)
 /*
	* Destructor, which cleans up after the interface Riemann solver class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CRoeRiemannSolver member functions.
//-----------------------------------------------------------------------------------


CRoeRiemannSolver::CRoeRiemannSolver
(
 CConfig *config_container
)
	:
		IRiemannSolver(config_container)
 /*
	* Constructor for a class based on Roe's Riemann solver.
	*/
{

}

//-----------------------------------------------------------------------------------

CRoeRiemannSolver::~CRoeRiemannSolver
(
 void
)
 /*
	* Destructor, which cleans up after Roe's Riemann solver class.
	*/
{

}


//-----------------------------------------------------------------------------------

void CRoeRiemannSolver::ComputeFlux
(
 const CMatrixAS3<as3double>     &wts,
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &solL,
 const CWorkMatrixAS3<as3double> &solR,
 CWorkMatrixAS3<as3double>       &flux
)
 /*
	* Function that cmoputes a unique flux state on a surface using Roe's method.
	*/
{
	// Loop over the integration points and compute the unique flux.
	for(size_t l=0; l<wts.size(); l++)
	{
		// The face is assumed to be taken from the perspective of the left solution.
		// Note, the residual is defined on the RHS: 
		//  dUDt = Res(vol) - Res(surf), thus the surf terms are -ve.
		const as3double wq = -C_HALF*wts[l]*met(0,l);
		const as3double nx =  met(1,l);
		const as3double ny =  met(2,l);

    // Compute the primitive data for the left state.
  	as3double tmp 	  	= C_ONE/solL(0,l);
  	const as3double	vxL = tmp*  solL(1,l);
  	const as3double vyL = tmp*  solL(2,l);
  	const as3double pL  = C_GM1*( solL(3,l) - C_HALF*( vxL*solL(1,l) + vyL*solL(2,l) ) );

  	// Compute the primitive data for the right state.
  	tmp                 = C_ONE/solR(0,l);
  	const as3double vxR = tmp*  solR(1,l);
  	const as3double vyR = tmp*  solR(2,l);
  	const as3double pR  = C_GM1*( solR(3,l) - C_HALF*( vxR*solR(1,l) + vyR*solR(2,l) ) );

  	// Compute the difference between the variables.
  	const as3double dr  = solR(0,l) - solL(0,l);
  	const as3double dru = solR(1,l) - solL(1,l);
  	const as3double drv = solR(2,l) - solL(2,l);
  	const as3double drE = solR(3,l) - solL(3,l);

  	// Compute the Roe average state.
  	const as3double zL = std::sqrt( solL(0,l) );
  	const as3double zR = std::sqrt( solR(0,l) );
  	tmp                = C_ONE/(zL + zR);

  	const as3double rHL = solL(3,l) + pL;
  	const as3double rHR = solR(3,l) + pR;

  	const as3double uAvg = tmp*(zL*vxL + zR*vxR);
  	const as3double vAvg = tmp*(zL*vyL + zR*vyR);
  	const as3double HAvg = tmp*(rHL/zL + rHR/zR);

    // Compute some abbreviations.
  	const as3double alphaAvg = C_HALF*(uAvg*uAvg + vAvg*vAvg);
  	tmp                      = C_GM1*(HAvg - alphaAvg);
  	const as3double a2Avg    = std::fabs(tmp);
  	const as3double aAvg     = std::sqrt(a2Avg);
  	const as3double vnAvg  	 = uAvg*nx + vAvg*ny;
  	const as3double ovaAvg   = C_ONE/aAvg;
  	const as3double ova2Avg  = C_ONE/a2Avg;

  	// Compute the absolute values of the eigenvalues.
  	as3double lam1 = std::fabs(vnAvg + aAvg);
  	as3double lam2 = std::fabs(vnAvg - aAvg);
  	as3double lam3 = std::fabs(vnAvg);

    // Apply the entropy fix.
    tmp  = mDelta*std::max(lam1, lam2);
  	lam1 = std::max(lam1, tmp);
  	lam2 = std::max(lam2, tmp);
  	lam3 = std::max(lam3, tmp);

  	// Some more abbreviations.
  	const as3double abv1 = C_HALF*(lam1 + lam2);
  	const as3double abv2 = C_HALF*(lam1 - lam2);
  	const as3double abv3 = abv1 - lam3;

  	const as3double abv4 = C_GM1*(alphaAvg*dr  - uAvg*dru - vAvg*drv + drE);
  	const as3double abv5 = nx*dru + ny*drv   - vnAvg*dr;
  	const as3double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
  	const as3double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;
	
  	const as3double vnL = vxL*nx + vyL*ny;
  	const as3double vnR = vxR*nx + vyR*ny;
  	const as3double pa  = pL + pR;

    // Compute the weighted Roe flux vector: 0.5*( FR + FL - |A|*(UR - UL) ).
		flux(0,l) = wq*(solL(0,l)*vnL + solR(0,l)*vnR
              - 	 (lam3*dr       + abv6));
    flux(1,l) = wq*(solL(1,l)*vnL + solR(1,l)*vnR + pa*nx
              -    (lam3*dru      + uAvg*abv6     + nx*abv7));
    flux(2,l) = wq*(solL(2,l)*vnL + solR(2,l)*vnR + pa*ny
              -    (lam3*drv      + vAvg*abv6     + ny*abv7));
    flux(3,l) = wq*(solL(3,l)*vnL + solR(3,l)*vnR + pL*vnL + pR*vnR
              -    (lam3*drE      + HAvg*abv6     + vnAvg*abv7));
	}
}







