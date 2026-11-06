#include "ic/CGaussianPressureIC.hpp"

#include "common/AS3Constants.hpp"
#include "CConfig.hpp"
#include "element/CPhysicalElement.hpp"
#include "geometry/CElementGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "solver/ISinglezoneSolver.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <ostream>


//-----------------------------------------------------------------------------------
// CGaussianPressureIC member functions.
//-----------------------------------------------------------------------------------


CGaussianPressureIC::CGaussianPressureIC
(
 const CConfig *config_container
)
	:
		IInitialCondition(config_container)
 /*
	* Constructor that reads the Gaussian pulse and freestream parameters.
	*/
{
	// Initialize the data, based on the user-specified values. Note, the unpacking
	// convention is based on the indices of mDataIC in IC_ExtractParam_GaussianPressure.
	mRatio    = config_container->GetDataIC(0);
	mWidth    = config_container->GetDataIC(1);
	mMachInf  = config_container->GetDataIC(2);
	mThetaInf = config_container->GetDataIC(3);
	mXCenter  = config_container->GetDataIC(4);
	mYCenter  = config_container->GetDataIC(5);

	// Default freestream pressure and temperature.
	mPressureInf    = static_cast<as3double>( 101325.0 );
	mTemperatureInf = static_cast<as3double>( 300.0    );

  // Deduce the freestream density.
  mDensityInf = mPressureInf/(C_RGAS*mTemperatureInf);

  // Compute the reference speed of sound.
  const as3double a    = std::sqrt(mPressureInf*C_GMA/mDensityInf);
  // Compute the flow speed.
  const as3double umag = mMachInf*a;
  
	// Compute the freestream velocity.
	mXVelocityInf = umag*std::cos(mThetaInf*C_PI/C_180);
  mYVelocityInf = umag*std::sin(mThetaInf*C_PI/C_180);
}

//-----------------------------------------------------------------------------------

CGaussianPressureIC::~CGaussianPressureIC
(
 void
)
 /*
	* Destructor for CGaussianPressureIC.
	*/
{

}

//-----------------------------------------------------------------------------------

void CGaussianPressureIC::InitializeSolution
(
 CConfig             *config_container,
 CSinglezoneGeometry *singlezone_geometry,
 ISinglezoneSolver  *solver_container
)
 /*
	* Function that initializes a Gaussian pressure pulse on the zone solution nodes.
	*/
{
	// Report output.
	std::cout << "    initial condition.... ";
	
	// Extract the total number of elements in this zone.
	const std::size_t nElem = singlezone_geometry->GetnElem();

	// Explicitly extract the required information.
	const as3double x0 = mXCenter;
	const as3double y0 = mYCenter;

	const as3double A0 = mRatio;
	const as3double b0 = mWidth;

	const as3double uInf   = mXVelocityInf;
	const as3double vInf   = mYVelocityInf;
	const as3double pInf   = mPressureInf;
	const as3double rhoInf = mDensityInf;

	// Abbreviations involving the width of the pressure.
	const as3double kappa  = std::log(C_TWO)/(b0*b0);

	// Abbreviation involving specific heat ratio.
	const as3double ovgm1  = C_ONE/C_GM1;

	// Loop over the elements and compute the solution.
	for(std::size_t i=0; i<nElem; i++)
	{
		// Extract solution in this element.
		auto& sol  = solver_container->GetPhysicalElement(i)->mSol2D;
		// Extract coordinates in this element.
		auto& coor = singlezone_geometry->GetElementGeometry(i)->GetCoordSolDOFs();

		// Loop over the internal solution nodes and compute the initial condition.
		for(std::size_t l=0; l<sol.col(); l++)
		{
			// Extract the coordinates.
			const as3double x = coor(0, l);
			const as3double y = coor(1, l);

			// Compute relative position w.r.t. pulse center.
			as3double dx = x - x0; 
			as3double dy = y - y0; 

			// Radial distance squared.
			const as3double r2 = dx*dx + dy*dy;

			// Compute the density, velocities and pressure.
			const as3double rho  = rhoInf;
			const as3double u    = uInf;
			const as3double v    = vInf;
			const as3double p    = pInf*( C_ONE + A0*std::exp(-kappa*r2) );
			// Compute the total energy.
			const as3double rhoE = p*ovgm1 + C_HALF*rho*( u*u + v*v );

			// Assemble conservative variables.
			sol(0, l) = rho;
			sol(1, l) = rho*u;
			sol(2, l) = rho*v;
			sol(3, l) = rhoE;
		}
	}
	// Report output.
	std::cout << "Done." << std::endl;
}
