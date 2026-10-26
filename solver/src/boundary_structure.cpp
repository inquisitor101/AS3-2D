#include "boundary_structure.hpp"

#include <cmath>
#include <cstddef>


//-----------------------------------------------------------------------------------
// IBoundary member functions.
//-----------------------------------------------------------------------------------


IBoundary::IBoundary
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *singlezone_geometry_container,
 const ISolver              *singlezone_solver_container,
 const CBoundaryFaceBatch   &boundary_batch,
 const IBoundaryParamMarker *param_bc
)
 /*
	* Constructor that selects the boundary Riemann solver and face kernels.
	*/
{
  // Extract the type of Riemann solver in each zone.
  const auto type_riemann = singlezone_solver_container->GetRiemannSolver()->GetTypeRiemannSolver(); 
 
  // Instantiate a Riemann solver on this interface.
	mRiemannSolverContainer = 
    CGenericFactory::CreateRiemannSolverContainer( config_container, type_riemann );

  // Initialize the compute kernels.
  InitializeComputeKernels( singlezone_solver_container, boundary_batch.GetFaceLocation() ); 
}

//-----------------------------------------------------------------------------------

IBoundary::~IBoundary
(
 void
)
 /*
	* Destructor for IBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void IBoundary::InitializeComputeKernels
(
 const ISolver *singlezone_solver_container,
 EFaceLocation  face_location
)
 /*
  * Function that selects the local face interpolation and residual kernels.
  */
{
  // Temporary lambda to determine which surface interpolation function to use.
  auto lGetFuncPointerInterpFace = [](auto* tensor_container, EFaceLocation face_location)
  {
	  // Create a function pointer for the iface in the imarker.
	  AInterpolateSurface FInterpFace;

	  // Definitions of the four interpolation functions on the (owner) iface, 
	  // which are given as lambda's that bind to std::function.
	  auto iSurfIMIN = 
      [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceIMIN(in...); };
	  auto iSurfIMAX = 
      [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceIMAX(in...); };
	  auto iSurfJMIN = 
      [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceJMIN(in...); };
	  auto iSurfJMAX = 
      [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceJMAX(in...); };

	  // Assign the appropriate interpolation functions for the (owner) iface.
	  switch(face_location)
	  {
	  	case(EFaceLocation::IMIN): {FInterpFace = iSurfIMIN; break;}
	  	case(EFaceLocation::IMAX): {FInterpFace = iSurfIMAX; break;}
	  	case(EFaceLocation::JMIN): {FInterpFace = iSurfJMIN; break;}
	  	case(EFaceLocation::JMAX): {FInterpFace = iSurfJMAX; break;}
	  	default: ERROR("Face is unknown.");
	  }

	  // Return the function pointer.
	  return FInterpFace;
  };


  // Temporary lambda to determine which surface residual function to use.
  auto lGetSurfaceResidualFace = [](auto* tensor_container, EFaceLocation face_location)
  {
	  // Create a function pointer for the iface in the imarker.
	  AComputeResidualFace FResFace;

	  // Definitions of the four interpolation functions on the (owner) iface, 
	  // which are given as lambda's that bind to std::function.
	  auto iSurfIMIN = 
      [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceIMIN(in...); };
	  auto iSurfIMAX = 
      [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceIMAX(in...); };
	  auto iSurfJMIN = 
      [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceJMIN(in...); };
	  auto iSurfJMAX = 
      [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceJMAX(in...); };

	  // Assign the appropriate interpolation functions for the (owner) iface.
	  switch(face_location)
	  {
	  	case(EFaceLocation::IMIN): {FResFace = iSurfIMIN; break;}
	  	case(EFaceLocation::IMAX): {FResFace = iSurfIMAX; break;}
	  	case(EFaceLocation::JMIN): {FResFace = iSurfJMIN; break;}
	  	case(EFaceLocation::JMAX): {FResFace = iSurfJMAX; break;}
	  	default: ERROR("Face is unknown.");
	  }

	  // Return the function pointer.
	  return FResFace;
  };

  // Obtain the relevant tensor-product container.
  auto* tensor_container = singlezone_solver_container->GetTensorProduct();

  // Assign the relevant interpolation function.
  mFInterpolateSurface = lGetFuncPointerInterpFace( tensor_container, face_location );

  // Assign the relevant residual function.
  mFComputeResidualFace = lGetSurfaceResidualFace( tensor_container, face_location );
}

//-----------------------------------------------------------------------------------

void IBoundary::ComputeBoundaryResidual
(
 const CSinglezoneGeometry *singlezone_geometry_container,
 ISolver                   *singlezone_solver_container,
 const CBoundaryFaceBatch   &face_batch,
 CFlattenedFaceIndex         face_info,
 CPoolMatrixAS3<as3double>  &workarray,
 as3double                   localtime
) const
 /*
  * Function that computes and lifts one physical boundary flux into its local face buffer.
  */
{
  // Extract the actual face.      
  const auto& boundary_face = face_batch.GetFace( face_info.mIndexFace );

  const auto iElem = boundary_face.mIndexElement;
  const EFaceLocation face_location = face_batch.GetFaceLocation();

  const auto* standard_element = singlezone_solver_container->GetStandardElement();

  // Get the respective number of integration points and weights on this face.
  const std::size_t nInt1D = standard_element->GetnInt1D();
  auto&        wInt1D = standard_element->GetwInt1D();
  // Also the number of solution variables.
  const std::size_t nVar   = singlezone_solver_container->GetnVar();

  // Get a pointer to the respective elements sharing this face.
  auto* elem = singlezone_solver_container->GetPhysicalElement(iElem);

  // Extract the metrics at the integration points on this face.
	auto& met = elem->GetSurfaceMetricInt(face_location);
	// Reference to the owner element solution.
	auto& sol = elem->mSol2D;	

  // Reference to the owner element residual.
  auto& res = elem->GetSurfaceResidual(face_location);
  res.reset();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> var_e = workarray.GetWorkMatrixAS3(nVar, nInt1D);
	CWorkMatrixAS3<as3double> var_b = workarray.GetWorkMatrixAS3(nVar, nInt1D);
	CWorkMatrixAS3<as3double> flux  = workarray.GetWorkMatrixAS3(nVar, nInt1D);


  // Compute the solution on the integration nodes of the owner element.
	mFInterpolateSurface(sol.data(), var_e.data(), nullptr, nullptr); 

  // Compute the boundary state, which is a slip-condition/symmetry/Euler wall.
  ComputeBoundaryState(met, var_e, var_b);

  // Compute the flux state, weighted by the integration nodes and metrics. 
  mRiemannSolverContainer->ComputeFlux(wInt1D, met, var_e, var_b, flux);

  // Compute the residual on the owned element.
  mFComputeResidualFace(flux.data(), nullptr, nullptr, res.data());
}


//-----------------------------------------------------------------------------------
// CSymmetryBoundary member functions.
//-----------------------------------------------------------------------------------


CSymmetryBoundary::CSymmetryBoundary
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *singlezone_geometry_container,
 const ISolver              *singlezone_solver_container,
 const CBoundaryFaceBatch   &boundary_batch,
 const IBoundaryParamMarker *param_bc
)
  : IBoundary(config_container,
              singlezone_geometry_container,
              singlezone_solver_container,
              boundary_batch,
              param_bc)
 /*
	* Constructor that initializes an inviscid symmetry boundary.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }
}

//-----------------------------------------------------------------------------------

CSymmetryBoundary::~CSymmetryBoundary
(
 void
)
 /*
	* Destructor for CSymmetryBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CSymmetryBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that reflects normal velocity to construct the exterior symmetry state.
  */
{
  for(std::size_t l=0; l<met.col(); l++)
  {
    const as3double nx = met(1,l);
		const as3double ny = met(2,l);

    const as3double ovrho = C_ONE / var_e(0,l);
          as3double u     = ovrho * var_e(1,l);
          as3double v     = ovrho * var_e(2,l);
    const as3double p     = C_GM1  * (   var_e(3,l) 
                          - C_HALF * ( u*var_e(1,l) + v*var_e(2,l) ) );

    // Compute the normal velocity component.
    const as3double un = u * nx + v * ny;

    // Subtract the normal components to obtain the tangential ones.
    u -= un * nx;
    v -= un * ny;
   
    // Assemble the boundary state.
    var_b(0,l) = var_e(0,l);
    var_b(1,l) = var_e(0,l) * u;
    var_b(2,l) = var_e(0,l) * v;
    var_b(3,l) = p * C_OVGM1 
               + C_HALF * var_e(0,l) * ( u*u + v*v );
  }
}


//-----------------------------------------------------------------------------------
// CRiemannStaticMassInletBoundary member functions.
//-----------------------------------------------------------------------------------


CRiemannStaticMassInletBoundary::CRiemannStaticMassInletBoundary
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *singlezone_geometry_container,
 const ISolver              *singlezone_solver_container,
 const CBoundaryFaceBatch   &boundary_batch,
 const IBoundaryParamMarker *param_bc
)
  : IBoundary(config_container,
              singlezone_geometry_container,
              singlezone_solver_container,
              boundary_batch,
              param_bc)
 /*
	* Constructor that stores the prescribed inlet density and velocity.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }

  // Convert after validating the boundary condition type.
  const auto* inlet_param =
    static_cast<const CStaticMassFlowInletRiemannParamMarker*>(param_bc);

  // Extract and assign the relevant prescribed values.
  mPrescribedDensity  = inlet_param->mDensity;
  mPrescribedVelocity = inlet_param->mVelocity;
}

//-----------------------------------------------------------------------------------

CRiemannStaticMassInletBoundary::~CRiemannStaticMassInletBoundary
(
 void
)
 /*
	* Destructor for CRiemannStaticMassInletBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CRiemannStaticMassInletBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that forms the exterior inlet state from prescribed density and velocity.
  */
{
  // Abbreviations.
  const as3double rb  = mPrescribedDensity;
  const as3double ub  = mPrescribedVelocity[0];
  const as3double vb  = mPrescribedVelocity[1];
  const as3double ekb = C_HALF * ( ub*ub + vb*vb );
  
  // Abbreviations involving gamma.
  const as3double gm1ovt  = C_HALF * C_GM1; 
  const as3double tgovgm1 = C_TWO * C_GMA * C_OVGM1;

  for(std::size_t l=0; l<met.col(); l++)
  {
    const as3double nx = met(1,l);
		const as3double ny = met(2,l);

    const as3double ovrho = C_ONE / var_e(0,l);
    const as3double u     = ovrho * var_e(1,l);
    const as3double v     = ovrho * var_e(2,l);
    const as3double p     = C_GM1  * (   var_e(3,l) 
                          - C_HALF * ( u*var_e(1,l) + v*var_e(2,l) ) );

    // Compute the acoustic information.
    const as3double a2 = C_GMA * p * ovrho;
    const as3double a  = std::sqrt( a2 );

    // Compute the normal velocity component.
    const as3double un = u * nx + v * ny;

    // Compute the normal velocity, based on the prescribed conditions.
    const as3double ubn = ub * nx + vb * ny; 
  
    // Compute the intermediate sound-speed ratio a_star/a along the interior
    // isentropic acoustic wave curve, using the prescribed normal velocity.
    const as3double ratio = C_ONE + gm1ovt * (un - ubn) / a;

    // Recover the corresponding static pressure. The boundary density is
    // prescribed independently when assembling the conservative state below.
    const as3double pb = p * std::pow(ratio, tgovgm1);

    // Assemble the boundary state.
    var_b(0,l) = rb;
    var_b(1,l) = rb * ub;
    var_b(2,l) = rb * vb;
    var_b(3,l) = pb * C_OVGM1 + rb * ekb;
  }
}


//-----------------------------------------------------------------------------------
// CRiemannPressureOutletBoundary member functions.
//-----------------------------------------------------------------------------------


CRiemannPressureOutletBoundary::CRiemannPressureOutletBoundary
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *singlezone_geometry_container,
 const ISolver              *singlezone_solver_container,
 const CBoundaryFaceBatch   &boundary_batch,
 const IBoundaryParamMarker *param_bc
)
  : IBoundary(config_container,
              singlezone_geometry_container,
              singlezone_solver_container,
              boundary_batch,
              param_bc)
 /*
	* Constructor that stores the prescribed outlet pressure.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }

  // Convert after validating the boundary condition type.
  const auto* outlet_param =
    static_cast<const CPressureOutletRiemannParamMarker*>(param_bc);

  // Extract and assign the relevant prescribed values.
  mPrescribedPressure = outlet_param->mPressure;
}

//-----------------------------------------------------------------------------------

CRiemannPressureOutletBoundary::~CRiemannPressureOutletBoundary
(
 void
)
 /*
	* Destructor for CRiemannPressureOutletBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CRiemannPressureOutletBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that forms the exterior outlet state with prescribed static pressure.
  */
{
  // Abbreviations.
  const as3double pb  = mPrescribedPressure;

  for(std::size_t l=0; l<met.col(); l++)
  {
    const as3double nx = met(1,l);
		const as3double ny = met(2,l);

    const as3double ovrho = C_ONE / var_e(0,l);
    const as3double u     = ovrho * var_e(1,l);
    const as3double v     = ovrho * var_e(2,l);
    const as3double p     = C_GM1  * (   var_e(3,l) 
                          - C_HALF * ( u*var_e(1,l) + v*var_e(2,l) ) );

    // Compute the acoustic information.
    const as3double a2 = C_GMA * p * ovrho;
    const as3double a  = std::sqrt( a2 );

    // Compute the normal velocity component.
    const as3double un = u * nx + v * ny;

    // Compute the Riemann invariant which is extrapolated.
    const as3double riemann = C_TWO * a * C_OVGM1 + un; 

    // Compute the specific entropy.
    const as3double s = p * std::pow( ovrho, C_GMA );

    // Compute the density, based on the extrapolated state.
    const as3double rb = std::pow( pb/s, C_OVGMA );

    // Compute the speed of sound, based on the extrapolated state.
    const as3double ab = std::sqrt( C_GMA * pb / rb );

    // Compute the normal velocity, based on the extrapolated state.
    const as3double ubn = riemann - C_TWO * ab * C_OVGM1; 
  
    // Compute the velocity of the extrapolated state.
    const as3double ub = u + (ubn - un) * nx;
    const as3double vb = v + (ubn - un) * ny;

    // Compute the specific kinetic energy of the extrapolated state.
    const as3double ekb = C_HALF * ( ub*ub + vb*vb );

    // Assemble the boundary state.
    var_b(0,l) = rb;
    var_b(1,l) = rb * ub;
    var_b(2,l) = rb * vb;
    var_b(3,l) = pb * C_OVGM1 + rb * ekb;
  }
}






