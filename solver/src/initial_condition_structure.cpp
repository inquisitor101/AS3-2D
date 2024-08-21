#include "initial_condition_structure.hpp"


//-----------------------------------------------------------------------------------
// IInitialCondition member functions.
//-----------------------------------------------------------------------------------

IInitialCondition::IInitialCondition
(
 CConfig *config_container
)
 /*
	* Constructor for the interface initial condition class.
	*/
{

}

//-----------------------------------------------------------------------------------

IInitialCondition::~IInitialCondition
(
 void
)
 /*
	* Destructor, which cleans up after the interface initial condition class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CGaussianPressureIC member functions.
//-----------------------------------------------------------------------------------


CGaussianPressureIC::CGaussianPressureIC
(
 CConfig *config_container
)
	:
		IInitialCondition(config_container)
 /*
	* Constructor for the Gaussian pressure initial condition class.
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
	mPressureInf    = 101325.0;
	mTemperatureInf = 300.0;

  // Deduce the freestream density.
  mDensityInf = mPressureInf/(GAS_CONSTANT*mTemperatureInf);

  // Compute the reference speed of sound.
  const as3double a    = std::sqrt(mPressureInf*GAMMA/mDensityInf);
  // Compute the flow speed.
  const as3double umag = mMachInf*a;
  
	// Compute the freestream velocity.
	mXVelocityInf = umag*std::cos(mThetaInf*PI_CONSTANT/180.0);
  mYVelocityInf = umag*std::sin(mThetaInf*PI_CONSTANT/180.0);
}

//-----------------------------------------------------------------------------------

CGaussianPressureIC::~CGaussianPressureIC
(
 void
)
 /*
	* Destructor, which cleans up after the gaussian pressure class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CGaussianPressureIC::InitializeSolution
(
 CConfig       *config_container,
 CZoneGeometry *zone_geometry,
 ISolver       *solver_container
)
 /*
	* Function that initializes a Gaussian pressure pulse in a zone.
	*/
{
	// Report output.
	std::cout << "    initial condition.... ";
	
	// Extract the total number of elements in this zone.
	const size_t nElem = zone_geometry->GetnElem();

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
	const as3double kappa  = std::log(2.0)/(b0*b0);

	// Abbreviation involving specific heat ratio.
	const as3double ovgm1  = 1.0/GAMMA_MINUS_ONE;

	// Loop over the elements and compute the solution.
	for(size_t i=0; i<nElem; i++)
	{
		// Extract solution in this element.
		auto& sol  = solver_container->GetPhysicalElement(i)->mSol2D;
		// Extract coordinates in this element.
		auto& coor = zone_geometry->GetElementGeometry(i)->GetCoordSolDOFs();

		// Loop over the internal solution nodes and compute the initial condition.
		for(size_t l=0; l<sol.col(); l++)
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
			const as3double p    = pInf*( 1.0 + A0*std::exp(-kappa*r2) );
			// Compute the total energy.
			const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

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


//-----------------------------------------------------------------------------------
// CIsentropicVortexIC member functions.
//-----------------------------------------------------------------------------------


CIsentropicVortexIC::CIsentropicVortexIC
(
 CConfig *config_container
)
	:
		IInitialCondition(config_container)
 /*
	* Constructor for the isentropic vortex initial condition class.
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
	mPressureInf    = 101325.0;
	mTemperatureInf = 300.0;

  // Deduce the freestream density.
  mDensityInf = mPressureInf/(GAS_CONSTANT*mTemperatureInf);

  // Compute the reference speed of sound.
  const as3double a    = std::sqrt(mPressureInf*GAMMA/mDensityInf);
  // Compute the flow speed.
  const as3double umag = mMachInf*a;
  
	// Compute the freestream velocity.
	mXVelocityInf = umag*std::cos(mThetaInf*PI_CONSTANT/180.0);
  mYVelocityInf = umag*std::sin(mThetaInf*PI_CONSTANT/180.0);
}

//-----------------------------------------------------------------------------------

CIsentropicVortexIC::~CIsentropicVortexIC
(
 void
)
 /*
	* Destructor, which cleans up after the isentropic vortex class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CIsentropicVortexIC::InitializeSolution
(
 CConfig       *config_container,
 CZoneGeometry *zone_geometry,
 ISolver       *solver_container
)
 /*
	* Function that initializes an isentropic vortex in a zone.
	*/
{
	// Report output.
	std::cout << "    initial condition.... ";
	
	// Extract the total number of elements in this zone.
	const size_t nElem = zone_geometry->GetnElem();

	// Explicitly extract the required information.
	const as3double x0 = mXCenter;
	const as3double y0 = mYCenter;

	const as3double A0 = mRatio;
	const as3double Rv = mWidth;

	const as3double uInf   = mXVelocityInf;
	const as3double vInf   = mYVelocityInf;
	const as3double pInf   = mPressureInf;
	const as3double rhoInf = mDensityInf;
	const as3double Tinf   = mTemperatureInf;
	const as3double aInf   = std::sqrt(pInf*GAMMA/rhoInf);
  
		// Abbreviations.
	const as3double ovgm1 =  1.0/GAMMA_MINUS_ONE;
  const as3double ovRv2 =  1.0/(Rv*Rv);
  const as3double alpha =  A0/(Rv*aInf);
  const as3double omega = -0.5*GAMMA*alpha*alpha;

	// Loop over the elements and compute the solution.
	for(size_t i=0; i<nElem; i++)
	{
		// Extract solution in this element.
		auto& sol  = solver_container->GetPhysicalElement(i)->mSol2D;
		// Extract coordinates in this element.
		auto& coor = zone_geometry->GetElementGeometry(i)->GetCoordSolDOFs();

		// Loop over the internal solution nodes and compute the initial condition.
		for(size_t l=0; l<sol.col(); l++)
		{
			// Extract the coordinates.
			const as3double x = coor(0, l);
			const as3double y = coor(1, l);

			// Compute relative position w.r.t. pulse center.
			as3double dx = x - x0; 
			as3double dy = y - y0; 

			// Radial distance squared.
			const as3double r2 = dx*dx + dy*dy;
 
			// Compute the stream function.
			const as3double psi = A0*std::exp(-0.5*r2*ovRv2);

			// Compute the velocities, pressure and density.
			const as3double u    = -ovRv2*dy*psi + uInf;
			const as3double v    =  ovRv2*dx*psi + vInf;
			const as3double p    =  pInf*std::exp( omega*std::exp(-r2*ovRv2) );
    	const as3double rho  =  p/(GAS_CONSTANT*Tinf);
			// Compute the total energy.
			const as3double rhoE = p*ovgm1 + 0.5*rho*( u*u + v*v );

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



