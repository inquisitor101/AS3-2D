#include "boundary/IBoundary.hpp"
#include "data/CPoolMatrixAS3.hpp"

#include <cstddef>

#include "compute_kernels/TensorProducts.hpp"
#include "CGenericFactory.hpp"
#include "geometry/face/CFlattenedFaceIndex.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "geometry/face/CBoundaryFaceBatch.hpp"
#include "element/CPhysicalElement.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "solver/ISinglezoneSolver.hpp"
#include "element/CStandardElement.hpp"


//-----------------------------------------------------------------------------------
// IBoundary member functions.
//-----------------------------------------------------------------------------------


IBoundary::IBoundary
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *singlezone_geometry_container,
 const ISinglezoneSolver    *singlezone_solver_container,
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
 const ISinglezoneSolver *singlezone_solver_container,
 EFaceLocation            face_location
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
 ISinglezoneSolver         *singlezone_solver_container,
 const CBoundaryFaceBatch  &face_batch,
 CFlattenedFaceIndex        face_info,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime
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
  auto&             wInt1D = standard_element->GetwInt1D();
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
