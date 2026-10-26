#include "interface_structure.hpp"

#include <cstddef>
#include <string>

#include <algorithm>


//-----------------------------------------------------------------------------------
// IInterface member functions.
//-----------------------------------------------------------------------------------


IInterface::IInterface
(
 const CConfig               *config_container,
 const CMultizoneSolver      *multizone_solver_container,
 const CInterfaceFaceBatch   &interface_batch
)
 /*
  * Constructor that records paired-zone metadata and creates interface operators.
  */
{
  // Assign the general interface properties.
  mIName = interface_batch.GetOwnerFaceName();
  mJName = interface_batch.GetNeighborFaceName();

  mIZone = interface_batch.GetOwnerZoneIndex();
  mJZone = interface_batch.GetNeighborZoneIndex();

  mIFace = interface_batch.GetOwnerFaceLocation();
  mJFace = interface_batch.GetNeighborFaceLocation();

  mNFace = interface_batch.GetnFace(); 

	mNPolyI = multizone_solver_container->GetSinglezoneSolver(mIZone)->GetStandardElement()->GetnPolySol();
	mNPolyJ = multizone_solver_container->GetSinglezoneSolver(mJZone)->GetStandardElement()->GetnPolySol();
}

//-----------------------------------------------------------------------------------

IInterface::~IInterface
(
 void
)
 /*
	* Destructor for IInterface.
	*/
{

}



//-----------------------------------------------------------------------------------
// CEEInterface member functions.
//-----------------------------------------------------------------------------------


CEEInterface::CEEInterface
(
 const CConfig               *config_container,
 const CMultizoneSolver      *multizone_solver_container,
 const CInterfaceFaceBatch   &interface_batch
)
  :
    IInterface(config_container,
               multizone_solver_container,
               interface_batch)
 /*
  * Constructor that selects Euler interface interpolation and residual kernels.
  */
{
  // Extract the relevant solvers.
  const auto* isolver = multizone_solver_container->GetSinglezoneSolver(mIZone);
  const auto* jsolver = multizone_solver_container->GetSinglezoneSolver(mJZone);

	// Extract the number of integration points in each marker zone.
	const std::size_t inint = isolver->GetStandardElement()->GetnInt1D();
	const std::size_t jnint = jsolver->GetStandardElement()->GetnInt1D();

	// Extract the type of DOFs in each zone.
	const ETypeDOF itype_dof = isolver->GetStandardElement()->GetTypeDOFsSol();
	const ETypeDOF jtype_dof = jsolver->GetStandardElement()->GetTypeDOFsSol();

	// For now, only use the same type of nodal points.
	if( itype_dof != jtype_dof )
  {
    ERROR("Currently, only same nodal points are supported.");
  }
  const ETypeDOF type_dof = itype_dof; // since i and j have the same distribution (for now).

	// Take the integration rule based on the highest polynomial.
	mNInt1D = std::max( inint, jnint );

	// Instantiate the appropriate (temporary) standard element containers in each zone.
	CStandardElement ielement(type_dof, mNPolyI, mNInt1D);
	CStandardElement jelement(type_dof, mNPolyJ, mNInt1D);

	// Obtain the integration weights on this interface.
	mWInt1D = ielement.GetwInt1D();

	// Instantiate the tensor-product containers in iZone and jZone.
	mITensorProductContainer = CGenericFactory::CreateTensorContainer( &ielement, isolver->GetnVar() );
	mJTensorProductContainer = CGenericFactory::CreateTensorContainer( &jelement, jsolver->GetnVar() );

	// Extract the type of Riemann solver in each zone.
	const auto iriemann = isolver->GetRiemannSolver()->GetTypeRiemannSolver(); 
	const auto jriemann = jsolver->GetRiemannSolver()->GetTypeRiemannSolver();
	
	// For now, force the Riemann solvers in both zones to be identical.
	if( iriemann != jriemann )
	{
		ERROR("For now, Riemann solvers at an interface must be identical.");
	}

	// Instantiate a Riemann solver on this interface.
	mRiemannSolverContainer = CGenericFactory::CreateRiemannSolverContainer( config_container, iriemann );

	// Ensure the number of working variables in the iZone is as expected.
	if( isolver->GetnVar() != jsolver->GetnVar() )
	{
		ERROR("Number of variables mismatches in iZone: " + std::to_string(mIZone) + " and jZone: " + std::to_string(mJZone) );
	}

  // Set the number of variables for this class.
  mNVar = isolver->GetnVar(); // since i and j have the same nVar (for now).


  // Initialize the compute kernerls for the solution interpolation and residual computation.
  InitializeComputeKernels();
}

//-----------------------------------------------------------------------------------

void CEEInterface::InitializeComputeKernels
(
 void
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
	  auto iSurfIMIN = [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceIMIN(in...); };
	  auto iSurfIMAX = [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceIMAX(in...); };
	  auto iSurfJMIN = [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceJMIN(in...); };
	  auto iSurfJMAX = [tensor_container](auto... in){ tensor_container->CompileTimeSurfaceJMAX(in...); };

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
	  auto iSurfIMIN = [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceIMIN(in...); };
	  auto iSurfIMAX = [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceIMAX(in...); };
	  auto iSurfJMIN = [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceJMIN(in...); };
	  auto iSurfJMAX = [tensor_container](auto... in){ tensor_container->CompileTimeResidualSurfaceJMAX(in...); };

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

  // Assign the relevant interpolation functions.
  mFInterpolateSurfaceI = lGetFuncPointerInterpFace( mITensorProductContainer.get(), mIFace );
  mFInterpolateSurfaceJ = lGetFuncPointerInterpFace( mJTensorProductContainer.get(), mJFace );

  // Assign the relevant residual functions.
  mFComputeResidualFaceI = lGetSurfaceResidualFace( mITensorProductContainer.get(), mIFace );
  mFComputeResidualFaceJ = lGetSurfaceResidualFace( mJTensorProductContainer.get(), mJFace );
}

//-----------------------------------------------------------------------------------

CEEInterface::~CEEInterface
(
 void
)
 /*
	* Destructor for CEEInterface.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEEInterface::ComputeInterfaceResidual
(
 CMultizoneSolver            *multizone_solver_container,
 const CInterfaceFaceBatch   &face_batch,
 CFlattenedFaceIndex          face_info,
 CPoolMatrixAS3<as3double>   &workarray,
 as3double                    localtime
) const
 /*
	* Function that computes and lifts one interface flux into both local face buffers.
	*/
{
  // Extract the relevant information.
  const std::size_t iZoneM = face_batch.GetOwnerZoneIndex();
  const std::size_t iZoneP = face_batch.GetNeighborZoneIndex();

  // Extract the actual face.
  const auto& interface_face = face_batch.GetFace( face_info.mIndexFace ); 

  const std::size_t iElemM = interface_face.mIndexElementOwner;
  const std::size_t iElemP = interface_face.mIndexElementNeighbor;

  const EFaceLocation face_location_m = face_batch.GetOwnerFaceLocation();
  const EFaceLocation face_location_p = face_batch.GetNeighborFaceLocation();

	// Get the solvers of this class.
	auto* solver_m = multizone_solver_container->GetSinglezoneSolver(iZoneM);
	auto* solver_p = multizone_solver_container->GetSinglezoneSolver(iZoneP);

  // Get the respective standard elements.
  const auto* standard_element_m = solver_m->GetStandardElement();
  const auto* standard_element_p = solver_p->GetStandardElement();

  // Get the respective number of integration points and weights on each face.
  const std::size_t nInt1D = standard_element_m->GetnInt1D();
  auto&        wInt1D = standard_element_m->GetwInt1D();
  // Also the number of solution variables.
  const std::size_t nVar   = solver_m->GetnVar();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> var_m = workarray.GetWorkMatrixAS3(nVar, nInt1D);
	CWorkMatrixAS3<as3double> var_p = workarray.GetWorkMatrixAS3(nVar, nInt1D);
	CWorkMatrixAS3<as3double> flux  = workarray.GetWorkMatrixAS3(nVar, nInt1D);

	// Get a pointer to the respective elements sharing this face.
	auto* elem_m = solver_m->GetPhysicalElement(iElemM);
	auto* elem_p = solver_p->GetPhysicalElement(iElemP);

	// Extract the metrics at the integration points on the owner face.
	auto& met_m = elem_m->GetSurfaceMetricInt(face_location_m);
	// Reference to the owner element solution.
	auto& sol_m = elem_m->mSol2D;	
  // Reference to the matching element solution.
	auto& sol_p = elem_p->mSol2D;
	
  // Reference to the owner element residual.
  auto& res_m = elem_m->GetSurfaceResidual(face_location_m);
  auto& res_p = elem_p->GetSurfaceResidual(face_location_p);
  res_m.reset();
  res_p.reset();
	

  // Compute the solution on the integration nodes of the owner element.
	mFInterpolateSurfaceI(sol_m.data(), var_m.data(), nullptr, nullptr); 

	// Compute the solution on the integration nodes of the matching element.
	mFInterpolateSurfaceJ(sol_p.data(), var_p.data(), nullptr, nullptr); 


  // Match the neighbor states to the owner's quadrature-point ordering.
  if( face_batch.IsNeighborReversed() )
  {
    for(std::size_t iVar=0; iVar<nVar; iVar++)
    {
      auto* row = var_p.data() + iVar*nInt1D;
      std::reverse(row, row+nInt1D);
    }
  }

  // Compute the flux state, weighted by the integration nodes and metrics. 
  mRiemannSolverContainer->ComputeFlux(wInt1D, met_m, var_m, var_p, flux);

  // Compute the residual on the owned element, which is on the iface boundary.
  mFComputeResidualFaceI(flux.data(), nullptr, nullptr, res_m.data());

	// For local conservation, negate the flux, since it leaves the owner element 
	// to enter the matching element.
	for(std::size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

  // Return the flux to the neighbor's native ordering before lifting it.
  if( face_batch.IsNeighborReversed() )
  {
    for(std::size_t iVar=0; iVar<nVar; iVar++)
    {
      auto* row = flux.data() + iVar*nInt1D;
      std::reverse(row, row+nInt1D);
    }
  }

	// Compute the residual on the matching element, which is on the jface boundary.
	mFComputeResidualFaceJ(flux.data(), nullptr, nullptr, res_p.data());
}





