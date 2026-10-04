#include "interface_structure.hpp"


//-----------------------------------------------------------------------------------
// IInterface member functions.
//-----------------------------------------------------------------------------------


IInterface::IInterface
(
 const CConfig               *config_container,
 const CMultizoneSolver      *multizone_solver_container,
 const CInterfaceFacesFamily &interface_family
)
 /*
  *
  */
{
  // Assign the general interface properties.
  mIName = interface_family.GetiFaceName();
  mJName = interface_family.GetjFaceName();

  mIZone = interface_family.GetiZone();
  mJZone = interface_family.GetjZone();

  mIFace = interface_family.GetiFaceLocation();
  mJFace = interface_family.GetjFaceLocation();

  mNFace = interface_family.GetnFace(); 

	mNPolyI = multizone_solver_container->GetSinglezoneSolver(mIZone)->GetStandardElement()->GetnPolySol();
	mNPolyJ = multizone_solver_container->GetSinglezoneSolver(mJZone)->GetStandardElement()->GetnPolySol();
}

//-----------------------------------------------------------------------------------

IInterface::~IInterface
(
 void
)
 /*
	* Destructor, which cleans up after the interface zone interface class.
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
 const CInterfaceFacesFamily &interface_family
)
  :
    IInterface(config_container,
               multizone_solver_container,
               interface_family)
 /*
  *
  */
{
  // Extract the relevant solvers.
  const auto* isolver = multizone_solver_container->GetSinglezoneSolver(mIZone);
  const auto* jsolver = multizone_solver_container->GetSinglezoneSolver(mJZone);

	// Extract the number of integration points in each marker zone.
	const unsigned short inint = isolver->GetStandardElement()->GetnInt1D();
	const unsigned short jnint = jsolver->GetStandardElement()->GetnInt1D();

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
  *
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
  mlInterpolateSurfaceI = lGetFuncPointerInterpFace( mITensorProductContainer.get(), mIFace );
  mlInterpolateSurfaceJ = lGetFuncPointerInterpFace( mJTensorProductContainer.get(), mJFace );

  // Assign the relevant residual functions.
  mlComputeResidualFaceI = lGetSurfaceResidualFace( mITensorProductContainer.get(), mIFace );
  mlComputeResidualFaceJ = lGetSurfaceResidualFace( mJTensorProductContainer.get(), mJFace );
}

//-----------------------------------------------------------------------------------

CEEInterface::~CEEInterface
(
 void
)
 /*
	* Destructor, which cleans up after the Euler equations zone interface class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEEInterface::ComputeInterfaceResidual
(
 CMultizoneSolver            *multizone_solver_container,
 const CInterfaceFacesFamily &family_face,
 CFlattenedFaceIndex          face_info,
 CPoolMatrixAS3<as3double>   &workarray,
 as3double                    localtime
) const
 /*
	* Function that computes the residual on a single element interface face. 
	*/
{
  // Extract the relevant information.
  const unsigned short iZoneM = family_face.GetiZone();
  const unsigned short iZoneP = family_face.GetjZone();

  // Extract the actual face.
  const auto& interface_face = family_face.GetInterfaceFace( face_info.mIndexFace ); 

  const size_t iElemM = interface_face.mIndexElementI;
  const size_t iElemP = interface_face.mIndexElementJ;

  const EFaceLocation face_location_m = family_face.GetiFaceLocation();
  const EFaceLocation face_location_p = family_face.GetjFaceLocation();

	// Get the solvers of this class.
	auto* solver_m = multizone_solver_container->GetSinglezoneSolver(iZoneM);
	auto* solver_p = multizone_solver_container->GetSinglezoneSolver(iZoneP);

  // Get the respective standard elements.
  const auto* standard_element_m = solver_m->GetStandardElement();
  const auto* standard_element_p = solver_p->GetStandardElement();

  // Get the respective number of integration points and weights on each face.
  const size_t nInt1D = standard_element_m->GetnInt1D();
  auto&        wInt1D = standard_element_m->GetwInt1D();
  // Also the number of solution variables.
  const size_t nVar   = solver_m->GetnVar();

// DEBUGGING start
  if( iZoneM != mIZone || iZoneP != mJZone )
  {
    std::cout << "iZoneM: " << iZoneM << ", mIZone: " << mIZone << ", "
              << "iZoneP: " << iZoneP << ", mJZone: " << mJZone << std::endl; 
  
    ERROR("zone indices do not match");
  }
  if( nInt1D != mNInt1D ) ERROR("nIntergration does not match.");
  if( face_location_m != mIFace || face_location_p != mJFace ) ERROR("face locations do not match.");
// DEBUGGING over

  // Consistency check.
  if( standard_element_m->GetnInt1D() != standard_element_p->GetnInt1D() 
      ||
      solver_m->GetnVar() != solver_p->GetnVar() )
  {
    ERROR("Interface contains faces with different properties.");
  }

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
	// Reference to the owner element residual.
  CMatrixAS3<as3double>* res_m;
  if( face_location_m == EFaceLocation::IMAX || face_location_m == EFaceLocation::JMAX )
  {
    res_m = &elem_m->mResMinus;
    res_m->reset();
  }
  else
  {
    res_m = &elem_m->mRes2D;
  }
  

	// Reference to the matching element solution.
	auto& sol_p = elem_p->mSol2D;
	// Reference to the matching element residual.
  CMatrixAS3<as3double>* res_p;
  if( face_location_p == EFaceLocation::IMAX || face_location_p == EFaceLocation::JMAX )
  {
    res_p = &elem_p->mResMinus;
    res_p->reset();
  }
  else
  {
    res_p = &elem_p->mRes2D;
  }

  // TODO: remove these 
  // Get relevant tensor products.
  auto* tensor_container_m = solver_m->GetTensorProduct();
  auto* tensor_container_p = solver_p->GetTensorProduct();

  // Get relevant Riemann solvers.
  auto* riemann_solver_m = solver_m->GetRiemannSolver();
  auto* riemann_solver_p = solver_p->GetRiemannSolver();


  // Consistency check.
  if( tensor_container_m->GetnVar() != tensor_container_p->GetnVar() )
  {
    ERROR("Number of variables in the tensor containers do not match.");
  }

  // Temporary lambda to determine which surface interpolation function to use.
  auto lGetFuncPointerInterpFace = [](auto* tensor_container, EFaceLocation face_location)
  {
	  // Create a function pointer for the iface in the imarker.
	  std::function<void(const as3double*, as3double*, as3double*, as3double*)> FInterpFace;

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

  auto InterpSurface_M = lGetFuncPointerInterpFace(tensor_container_m, face_location_m);
  auto InterpSurface_P = lGetFuncPointerInterpFace(tensor_container_p, face_location_p);


  // Temporary lambda to determine which surface residual function to use.
  auto lGetSurfaceResidualFace = [](auto* tensor_container, EFaceLocation face_location)
  {
	  // Create a function pointer for the iface in the imarker.
	  std::function<void(const as3double*, as3double*, as3double*, as3double*)> FResFace;

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

  auto ComputeResFace_M = lGetSurfaceResidualFace(tensor_container_m, face_location_m);
  auto ComputeResFace_P = lGetSurfaceResidualFace(tensor_container_p, face_location_p);


	// Compute the solution on the integration nodes of the owner element.
	InterpSurface_M(sol_m.data(), var_m.data(), nullptr, nullptr); 

	// Compute the solution on the integration nodes of the matching element.
	InterpSurface_P(sol_p.data(), var_p.data(), nullptr, nullptr); 



  // Compute the flux state, weighted by the integration nodes and metrics. 
  riemann_solver_m->ComputeFlux(wInt1D, met_m, var_m, var_p, flux);

  // Compute the residual on the owned element, which is on the iface boundary.
  ComputeResFace_M(flux.data(), nullptr, nullptr, res_m->data());

	// For local conservation, negate the flux, since it leaves the owner element 
	// to enter the matching element.
	for(size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

	// Compute the residual on the matching element, which is on the jface boundary.
	ComputeResFace_P(flux.data(), nullptr, nullptr, res_p->data());


  // CHANGED: new version is below..

	//// Compute the solution on the integration nodes of the owner element.
	//mlInterpolateSurfaceI(sol_m.data(), var_m.data(), nullptr, nullptr); 

	//// Compute the solution on the integration nodes of the matching element.
	//mlInterpolateSurfaceJ(sol_p.data(), var_p.data(), nullptr, nullptr); 



  //// Compute the flux state, weighted by the integration nodes and metrics. 
  //mRiemannSolverContainer->ComputeFlux(wInt1D, met_m, var_m, var_p, flux);

  //// Compute the residual on the owned element, which is on the iface boundary.
  //mlComputeResidualFaceI(flux.data(), nullptr, nullptr, res_m->data());

	//// For local conservation, negate the flux, since it leaves the owner element 
	//// to enter the matching element.
	//for(size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

	//// Compute the residual on the matching element, which is on the jface boundary.
	//mlComputeResidualFaceJ(flux.data(), nullptr, nullptr, res_p->data());
}





