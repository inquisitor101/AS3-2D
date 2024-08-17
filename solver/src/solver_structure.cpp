#include "solver_structure.hpp"


//-----------------------------------------------------------------------------------
// ISolver member functions.
//-----------------------------------------------------------------------------------


ISolver::ISolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
	:
		mZoneID(iZone)
 /*
	* Constructor for the interface solver class.
	*/
{
	// Instantiate a standard element object.
	mStandardElementContainer = CGenericFactory::CreateStandardElement(config_container, iZone);

	// Instantiate a tensor-product object.
	mTensorProductContainer = CGenericFactory::CreateTensorContainer(mStandardElementContainer.get());

	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(iZone);

	// Ensure the zone matches this one.
	if( iZone != zone->GetZoneID() ) ERROR("Geometry zone does not match Solver zone.");
}

//-----------------------------------------------------------------------------------

ISolver::~ISolver
(
 void
)
 /*
	* Destructor, which cleans up after the interface solver class.
	*/
{

}


//-----------------------------------------------------------------------------------
// CEESolver member functions.
//-----------------------------------------------------------------------------------


CEESolver::CEESolver
(
 CConfig       *config_container,
 CGeometry     *geometry_container,
 unsigned short iZone
)
	:
		ISolver(config_container, geometry_container, iZone)
 /*
	* Constructor for the (non-linear) Euler equations class.
	*/
{

}

//-----------------------------------------------------------------------------------

CEESolver::~CEESolver
(
 void
)
 /*
	* Destructor, which cleans up after the Euler equations class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEESolver::InitPhysicalElements
(
 CConfig   *config_container,
 CGeometry *geometry_container
)
 /*
	* Function that initializes the physical elements. 
	*/
{
	// Get a reference to the current zone.
	auto* zone = geometry_container->GetZoneGeometry(mZoneID);

	// Allocate memory for the physical elements in this solver.
	mPhysicalElementContainer.resize( zone->GetnElem() );
	
	// Instantiate physical element objects. The reason this needs a factory is because we might have 
	// a specialized solver (e.g. EE), but with a buffer zone, such as a PML. Hence, the elements
	// require additional variables. For now, avoid an interface class, as this object is the most 
	// performance-reliant data class in the code. Perhaps it is better to keep a single class, but 
	// instantiate it differently, depending on the type of solver and buffer layer specification.
	for(size_t i=0; i<mPhysicalElementContainer.size(); i++)
	{
		mPhysicalElementContainer[i] = CGenericFactory::CreatePhysicalElement(config_container,
				                                                                  mStandardElementContainer.get(),
																																					mTensorProductContainer.get(),
				                                                                  zone->GetElementGeometry(i),
																																					mZoneID, mNVar);
	}
}

//-----------------------------------------------------------------------------------

void CEESolver::ComputeVolumeResidual
(
 size_t                     iElem,
 as3double                  localtime,
 as3vector1d<as3double>    &monitordata,
 CPoolMatrixAS3<as3double> &workarray
)
 /*
	* Function that computes the volume terms in a EE-type PDE. 
	*/
{
	// Extract the number of integration points in 2D.
	size_t nInt2D = mStandardElementContainer->GetnInt2D();
	// Extract the metrics at the volume solution points.
	auto& metrics = mPhysicalElementContainer[iElem]->mMetricInt2D;
	// Extract the quadrature integration weights on the standard element.
	auto& wInt2D  = mStandardElementContainer->GetwInt2D();
	// Reference to the current element solution.
	auto& sol     = mPhysicalElementContainer[iElem]->mSol2D;
	// Reference to the current element residual.
	auto& res     = mPhysicalElementContainer[iElem]->mRes2D;


	// Borrow memory for the solution and its gradient.
	CWorkMatrixAS3<as3double> Var    = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
	CWorkMatrixAS3<as3double> dVarDx = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
	CWorkMatrixAS3<as3double> dVarDy = workarray.GetWorkMatrixAS3(mNVar, nInt2D);


	// Compute the solution and its (parametric) gradient at the volume integration points.
	mTensorProductContainer->Volume(mNVar, sol.data(),
			                            Var.data(), dVarDx.data(), dVarDy.data());

	// Convert the gradient from parametric to Cartesian coordinates.
	mPhysicalElementContainer[iElem]->ConvertGradParamToCartVolInt(dVarDx, dVarDy);

	// Abbreviations.
	const as3double gm1 = GAMMA_MINUS_ONE;

	// Loop over the integration points and compute the weighted terms.
	for(size_t l=0; l<nInt2D; l++)
	{
		// Compute the weighting factor: Jacobian multiplied by the integration weight.
		const as3double weight = metrics(0,l)*wInt2D[l];

		// Assemble the relevant (weighted) metrics in the x- and y-derivatives.
		const as3double wdrdx = weight*metrics(1,l);
		const as3double wdrdy = weight*metrics(2,l);
		const as3double wdsdx = weight*metrics(3,l);
		const as3double wdsdy = weight*metrics(4,l);

    // Compute the primitive variables.
    const as3double rho   = Var(0,l);
    const as3double ovrho = 1.0/rho;
    const as3double u     = ovrho*Var(1,l);
    const as3double v     = ovrho*Var(2,l);
    const as3double p     = gm1*( Var(3,l) - 0.5*(u*Var(1,l) + v*Var(2,l)) );

    // Compute the inviscid flux in the x-direction. 
    const as3double fx0 =     Var(1,l);
    const as3double fx1 = u*  Var(1,l) + p;
    const as3double fx2 = v*  Var(1,l);
    const as3double fx3 = u*( Var(3,l) + p );

    // Compute the inviscid flux in the y-direction. 
    const as3double fy0 =     Var(2,l);
    const as3double fy1 = u*  Var(2,l);
    const as3double fy2 = v*  Var(2,l) + p;
    const as3double fy3 = v*( Var(3,l) + p );

		// Compute the residual flux terms in the r-parametric direction.
		dVarDx(0,l) = fx0*wdrdx + fy0*wdrdy;
		dVarDx(1,l) = fx1*wdrdx + fy1*wdrdy;
		dVarDx(2,l) = fx2*wdrdx + fy2*wdrdy;
		dVarDx(3,l) = fx3*wdrdx + fy3*wdrdy;

		// Compute the residual flux terms in the s-parametric direction.
		dVarDy(0,l) = fx0*wdsdx + fy0*wdsdy;
		dVarDy(1,l) = fx1*wdsdx + fy1*wdsdy;
		dVarDy(2,l) = fx2*wdsdx + fy2*wdsdy;
		dVarDy(3,l) = fx3*wdsdx + fy3*wdsdy;
	}


	// Scatter the volume residual terms back to the actual residual. Note, this resets the residual! 
	mTensorProductContainer->ResidualVolume(mNVar, 
			                                    nullptr, 
																					dVarDx.data(), 
																					dVarDy.data(), 
																					res.data());
}







