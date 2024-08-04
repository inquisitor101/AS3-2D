#include "physical_element_structure.hpp"


//-----------------------------------------------------------------------------------
// CPhysicalElement member functions.
//-----------------------------------------------------------------------------------


CPhysicalElement::CPhysicalElement
(
 CConfig          *config_container,
 CStandardElement *standard_element,
 ITensorProduct   *tensor_container,
 CElementGeometry *element_geometry,
 unsigned short    iZone,
 unsigned short    nVar
)
	:
		mNPoly( config_container->GetnPolySol(iZone) )
 /*
	* Constructor for the physical element, which initiates an element in physical space.
	*/
{
	// Deduce the number of points in 1D and 2D.
	const size_t nSol1D = mNPoly+1;
	const size_t nSol2D = nSol1D*nSol1D;

	// Allocate memory for the solution in 2D.
	mSol2D.resize( nVar, nSol2D );

	
	// Compute the metrics at the volume solution points.
	ComputeMetricsSolVolume(standard_element, tensor_container, element_geometry->GetCoordSolDOFs());

	// Compute the metrics at the volume integration points.
	ComputeMetricsIntVolume(standard_element, tensor_container, element_geometry->GetCoordSolDOFs());

	// Compute the metrics at the integration points of the four surfaces of the quadrilateral element.
	ComputeMetricsIntSurfIMIN(standard_element, tensor_container, element_geometry->GetCoordSolDOFs());
	ComputeMetricsIntSurfIMAX(standard_element, tensor_container, element_geometry->GetCoordSolDOFs());
	ComputeMetricsIntSurfJMIN(standard_element, tensor_container, element_geometry->GetCoordSolDOFs());
	ComputeMetricsIntSurfJMAX(standard_element, tensor_container, element_geometry->GetCoordSolDOFs());
}

//-----------------------------------------------------------------------------------

CPhysicalElement::~CPhysicalElement
(
 void
)
 /*
	* Destructor, which cleans up after the physical element class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CPhysicalElement::ComputeMetricsSolVolume
(
 CStandardElement      *standard_element,
 ITensorProduct        *tensor_container, 
 CMatrixAS3<as3double> &coord
)
 /*
	* Function that computes the metrics at the volume solution points.
	*/
{
	// Dimension of coordinates.
	const size_t nDim   = 2;
	// Number of items in the metrics. See header for info on those items.
	const size_t nItem  = 5;
	// Number of volume solution points in 1D.
	const size_t nSol1D = standard_element->GetnSol1D(); 
	// Number of volume solution points in 2D.
	const size_t nSol2D = standard_element->GetnSol2D(); 

	// Define constants, such that precision is implicit.
	const as3double one  = (as3double) 1.0;
	const as3double zero = (as3double) 0.0;

	// Allocate memory for the volume metrics.
	mMetricSol2D.resize( nItem, nSol2D );

	// Allocate memory for the derivative of the solution at the solution points.
	CMatrixAS3<as3double> dSolDr( nDim, nSol2D );
	CMatrixAS3<as3double> dSolDs( nDim, nSol2D );

	// Create an identity matrix that serves as the Lagrange interpolation on the same solution points.
	CMatrixAS3<as3double> identity = NLinearAlgebra::CreateIdentityMatrix<as3double>(nSol1D);

	// Compute the derivatives on the standard element.
	tensor_container->CustomVolume( nSol1D, nDim, nSol1D,
			                            identity.data(), standard_element->GetDerLagrangeSol1DTrans().data(),
			                            coord.data(), nullptr,
																	dSolDr.data(), dSolDs.data());


	// Compute the metrics.
	for(size_t i=0; i<nSol2D; i++)
	{
		const as3double dxdr = dSolDr(0, i);
		const as3double dxds = dSolDs(0, i);
		const as3double dydr = dSolDr(1, i);
		const as3double dyds = dSolDs(1, i);

		// Jacobian of the transformation.
		const as3double J    = dxdr*dyds - dydr*dxds;

		// Check if the element makes sense.
		if( J <= zero ) ERROR("Negative Jacobian is encountered.");

		// Inverse of the Jacobian.
		const as3double Jinv = one/J;

		// Compute and store the metrics.
		mMetricSol2D(0, i) =  J;             // |J|
		mMetricSol2D(1, i) =  Jinv*( dyds ); // drdx
		mMetricSol2D(2, i) = -Jinv*( dxds ); // drdy
		mMetricSol2D(3, i) = -Jinv*( dydr ); // dsdx
		mMetricSol2D(4, i) =  Jinv*( dxdr ); // dsdy
	}
}

//-----------------------------------------------------------------------------------

void CPhysicalElement::ComputeMetricsIntVolume
(
 CStandardElement      *standard_element,
 ITensorProduct        *tensor_container, 
 CMatrixAS3<as3double> &coord
)
 /*
	* Function that computes the metrics at the volume integration points.
	*/
{
	// Dimension of coordinates.
	const size_t nDim   = 2;
	// Number of items in the metrics. See header for info on those items.
	const size_t nItem  = 5;
	// Number of volume integration points in 2D.
	const size_t nInt2D = standard_element->GetnInt2D(); 

	// Define constants, such that precision is implicit.
	const as3double one  = (as3double) 1.0;
	const as3double zero = (as3double) 0.0;

	// Allocate memory for the volume metrics.
	mMetricInt2D.resize( nItem, nInt2D );

	// Allocate memory for the derivative of the solution at the integration points.
	CMatrixAS3<as3double> dSolDr( nDim, nInt2D );
	CMatrixAS3<as3double> dSolDs( nDim, nInt2D );

	// Compute the derivatives on the standard element.
	tensor_container->Volume(nDim, coord.data(), nullptr,
			                     dSolDr.data(), dSolDs.data());

	// Compute the metrics.
	for(size_t i=0; i<nInt2D; i++)
	{
		const as3double dxdr = dSolDr(0, i);
		const as3double dxds = dSolDs(0, i);
		const as3double dydr = dSolDr(1, i);
		const as3double dyds = dSolDs(1, i);

		// Jacobian of the transformation.
		const as3double J    = dxdr*dyds - dydr*dxds;

		// Check if the element makes sense.
		if( J <= zero ) ERROR("Negative Jacobian is encountered.");

		// Inverse of the Jacobian.
		const as3double Jinv = one/J;

		// Compute and store the metrics.
		mMetricInt2D(0, i) =  J;             // |J|
		mMetricInt2D(1, i) =  Jinv*( dyds ); // drdx
		mMetricInt2D(2, i) = -Jinv*( dxds ); // drdy
		mMetricInt2D(3, i) = -Jinv*( dydr ); // dsdx
		mMetricInt2D(4, i) =  Jinv*( dxdr ); // dsdy
	}
}

//-----------------------------------------------------------------------------------

void CPhysicalElement::ComputeMetricsIntSurfIMIN
(
 CStandardElement      *standard_element,
 ITensorProduct        *tensor_container, 
 CMatrixAS3<as3double> &coord
)
 /*
	* Function that computes the metrics at the IMIN surface integration points.
	*/
{
	// Dimension of coordinates.
	const size_t nDim   = 2;
	// Number of items in the metrics. See header for info on those items.
	const size_t nItem  = 7;
	// Number of volume integration points in 1D.
	const size_t nInt1D = standard_element->GetnInt1D(); 

	// Define constants, such that precision is implicit.
	const as3double one  = (as3double) 1.0;
	const as3double zero = (as3double) 0.0;

	// Allocate memory for the surface metrics.
	mMetricIntIMin1D.resize( nItem, nInt1D );

	// Allocate memory for the derivative of the solution at the integration points.
	CMatrixAS3<as3double> dSolDr( nDim, nInt1D );
	CMatrixAS3<as3double> dSolDs( nDim, nInt1D );

	// Compute the derivatives on the standard element.
	tensor_container->SurfaceIMIN(nDim, coord.data(), nullptr,
			                     dSolDr.data(), dSolDs.data());

	// Compute the metrics.
	for(size_t i=0; i<nInt1D; i++)
	{
		const as3double dxdr = dSolDr(0, i);
		const as3double dxds = dSolDs(0, i);
		const as3double dydr = dSolDr(1, i);
		const as3double dyds = dSolDs(1, i);

		// Jacobian of the transformation.
		const as3double J    = dxdr*dyds - dydr*dxds;

		// Check if the element makes sense.
		if( J <= zero ) ERROR("Negative Jacobian is encountered.");

		// Inverse of the Jacobian.
		const as3double Jinv = one/J;

		// Compute the normal in the negative r-direction. Note, IMIN surface normal points inward.
		const as3double nx = -dyds;
		const as3double ny =  dxds;
	
		// Compute the magnitude of this normal.
		const as3double ll = std::sqrt( nx*nx + ny*ny );

		// Compute and store the metrics.
		mMetricIntIMin1D(0, i) =  ll;            // norm(n)
		mMetricIntIMin1D(1, i) =  nx/ll;         // nx/norm(n)
		mMetricIntIMin1D(2, i) =  ny/ll;         // ny/norm(n)

		mMetricIntIMin1D(3, i) =  Jinv*( dyds ); // drdx
		mMetricIntIMin1D(4, i) = -Jinv*( dxds ); // drdy
		mMetricIntIMin1D(5, i) = -Jinv*( dydr ); // dsdx
		mMetricIntIMin1D(6, i) =  Jinv*( dxdr ); // dsdy
	}
}

//-----------------------------------------------------------------------------------

void CPhysicalElement::ComputeMetricsIntSurfIMAX
(
 CStandardElement      *standard_element,
 ITensorProduct        *tensor_container, 
 CMatrixAS3<as3double> &coord
)
 /*
	* Function that computes the metrics at the IMAX surface integration points.
	*/
{
	// Dimension of coordinates.
	const size_t nDim   = 2;
	// Number of items in the metrics. See header for info on those items.
	const size_t nItem  = 7;
	// Number of volume integration points in 1D.
	const size_t nInt1D = standard_element->GetnInt1D(); 

	// Define constants, such that precision is implicit.
	const as3double one  = (as3double) 1.0;
	const as3double zero = (as3double) 0.0;

	// Allocate memory for the surface metrics.
	mMetricIntIMax1D.resize( nItem, nInt1D );

	// Allocate memory for the derivative of the solution at the integration points.
	CMatrixAS3<as3double> dSolDr( nDim, nInt1D );
	CMatrixAS3<as3double> dSolDs( nDim, nInt1D );

	// Compute the derivatives on the standard element.
	tensor_container->SurfaceIMAX(nDim, coord.data(), nullptr,
			                     dSolDr.data(), dSolDs.data());

	// Compute the metrics.
	for(size_t i=0; i<nInt1D; i++)
	{
		const as3double dxdr = dSolDr(0, i);
		const as3double dxds = dSolDs(0, i);
		const as3double dydr = dSolDr(1, i);
		const as3double dyds = dSolDs(1, i);

		// Jacobian of the transformation.
		const as3double J    = dxdr*dyds - dydr*dxds;

		// Check if the element makes sense.
		if( J <= zero ) ERROR("Negative Jacobian is encountered.");

		// Inverse of the Jacobian.
		const as3double Jinv = one/J;

		// Compute the normal in the positive r-direction. Note, IMAX surface normal points outward.
		const as3double nx =  dyds;
		const as3double ny = -dxds;
	
		// Compute the magnitude of this normal.
		const as3double ll = std::sqrt( nx*nx + ny*ny );

		// Compute and store the metrics.
		mMetricIntIMax1D(0, i) =  ll;            // norm(n)
		mMetricIntIMax1D(1, i) =  nx/ll;         // nx/norm(n)
		mMetricIntIMax1D(2, i) =  ny/ll;         // ny/norm(n)

		mMetricIntIMax1D(3, i) =  Jinv*( dyds ); // drdx
		mMetricIntIMax1D(4, i) = -Jinv*( dxds ); // drdy
		mMetricIntIMax1D(5, i) = -Jinv*( dydr ); // dsdx
		mMetricIntIMax1D(6, i) =  Jinv*( dxdr ); // dsdy
	}
}

//-----------------------------------------------------------------------------------

void CPhysicalElement::ComputeMetricsIntSurfJMIN
(
 CStandardElement      *standard_element,
 ITensorProduct        *tensor_container, 
 CMatrixAS3<as3double> &coord
)
 /*
	* Function that computes the metrics at the JMIN surface integration points.
	*/
{
	// Dimension of coordinates.
	const size_t nDim   = 2;
	// Number of items in the metrics. See header for info on those items.
	const size_t nItem  = 7;
	// Number of volume integration points in 1D.
	const size_t nInt1D = standard_element->GetnInt1D(); 

	// Define constants, such that precision is implicit.
	const as3double one  = (as3double) 1.0;
	const as3double zero = (as3double) 0.0;

	// Allocate memory for the surface metrics.
	mMetricIntJMin1D.resize( nItem, nInt1D );

	// Allocate memory for the derivative of the solution at the integration points.
	CMatrixAS3<as3double> dSolDr( nDim, nInt1D );
	CMatrixAS3<as3double> dSolDs( nDim, nInt1D );

	// Compute the derivatives on the standard element.
	tensor_container->SurfaceJMIN(nDim, coord.data(), nullptr,
			                     dSolDr.data(), dSolDs.data());

	// Compute the metrics.
	for(size_t i=0; i<nInt1D; i++)
	{
		const as3double dxdr = dSolDr(0, i);
		const as3double dxds = dSolDs(0, i);
		const as3double dydr = dSolDr(1, i);
		const as3double dyds = dSolDs(1, i);

		// Jacobian of the transformation.
		const as3double J    = dxdr*dyds - dydr*dxds;

		// Check if the element makes sense.
		if( J <= zero ) ERROR("Negative Jacobian is encountered.");

		// Inverse of the Jacobian.
		const as3double Jinv = one/J;

		// Compute the normal in the negative s-direction. Note, JMIN surface normal points inward.
		const as3double nx =  dydr;
		const as3double ny = -dxdr;
	
		// Compute the magnitude of this normal.
		const as3double ll = std::sqrt( nx*nx + ny*ny );

		// Compute and store the metrics.
		mMetricIntJMin1D(0, i) =  ll;            // norm(n)
		mMetricIntJMin1D(1, i) =  nx/ll;         // nx/norm(n)
		mMetricIntJMin1D(2, i) =  ny/ll;         // ny/norm(n)

		mMetricIntJMin1D(3, i) =  Jinv*( dyds ); // drdx
		mMetricIntJMin1D(4, i) = -Jinv*( dxds ); // drdy
		mMetricIntJMin1D(5, i) = -Jinv*( dydr ); // dsdx
		mMetricIntJMin1D(6, i) =  Jinv*( dxdr ); // dsdy
	}
}

//-----------------------------------------------------------------------------------

void CPhysicalElement::ComputeMetricsIntSurfJMAX
(
 CStandardElement      *standard_element,
 ITensorProduct        *tensor_container, 
 CMatrixAS3<as3double> &coord
)
 /*
	* Function that computes the metrics at the JMAX surface integration points.
	*/
{
	// Dimension of coordinates.
	const size_t nDim   = 2;
	// Number of items in the metrics. See header for info on those items.
	const size_t nItem  = 7;
	// Number of volume integration points in 1D.
	const size_t nInt1D = standard_element->GetnInt1D(); 

	// Define constants, such that precision is implicit.
	const as3double one  = (as3double) 1.0;
	const as3double zero = (as3double) 0.0;

	// Allocate memory for the surface metrics.
	mMetricIntJMax1D.resize( nItem, nInt1D );

	// Allocate memory for the derivative of the solution at the integration points.
	CMatrixAS3<as3double> dSolDr( nDim, nInt1D );
	CMatrixAS3<as3double> dSolDs( nDim, nInt1D );

	// Compute the derivatives on the standard element.
	tensor_container->SurfaceJMAX(nDim, coord.data(), nullptr,
			                     dSolDr.data(), dSolDs.data());

	// Compute the metrics.
	for(size_t i=0; i<nInt1D; i++)
	{
		const as3double dxdr = dSolDr(0, i);
		const as3double dxds = dSolDs(0, i);
		const as3double dydr = dSolDr(1, i);
		const as3double dyds = dSolDs(1, i);

		// Jacobian of the transformation.
		const as3double J    = dxdr*dyds - dydr*dxds;

		// Check if the element makes sense.
		if( J <= zero ) ERROR("Negative Jacobian is encountered.");

		// Inverse of the Jacobian.
		const as3double Jinv = one/J;

		// Compute the normal in the positive s-direction. Note, JMAX surface normal points outward.
		const as3double nx = -dydr;
		const as3double ny =  dxdr;
	
		// Compute the magnitude of this normal.
		const as3double ll = std::sqrt( nx*nx + ny*ny );

		// Compute and store the metrics.
		mMetricIntJMax1D(0, i) =  ll;            // norm(n)
		mMetricIntJMax1D(1, i) =  nx/ll;         // nx/norm(n)
		mMetricIntJMax1D(2, i) =  ny/ll;         // ny/norm(n)

		mMetricIntJMax1D(3, i) =  Jinv*( dyds ); // drdx
		mMetricIntJMax1D(4, i) = -Jinv*( dxds ); // drdy
		mMetricIntJMax1D(5, i) = -Jinv*( dydr ); // dsdx
		mMetricIntJMax1D(6, i) =  Jinv*( dxdr ); // dsdy
	}
}




