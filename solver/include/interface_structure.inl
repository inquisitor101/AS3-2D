

//-----------------------------------------------------------------------------------
// Implementation of the inlined functions in IInterface. 
//-----------------------------------------------------------------------------------


auto IInterface::GetFuncPointerInterpFaceI
(
 void
)
 /*
	* Function which returns a function pointer for the interpolation on 
	* the (owner) face belonging to the i-marker.
	*/
{
	// Create a function pointer for the iface in the imarker.
	std::function<void(const size_t, const as3double*, as3double*, as3double*, as3double*)> FInterpFaceI;

	// Definitions of the four interpolation functions on the (owner) iface, 
	// which are given as lambda's that bind to std::function.
	auto iSurfIMIN = [this](auto... in){ mITensorProductContainer->SurfaceIMIN(in...); };
	auto iSurfIMAX = [this](auto... in){ mITensorProductContainer->SurfaceIMAX(in...); };
	auto iSurfJMIN = [this](auto... in){ mITensorProductContainer->SurfaceJMIN(in...); };
	auto iSurfJMAX = [this](auto... in){ mITensorProductContainer->SurfaceJMAX(in...); };

	// Assign the appropriate interpolation functions for the (owner) iface.
	switch(mIFace)
	{
		case(EFaceLocation::IMIN): {FInterpFaceI = iSurfIMIN; break;}
		case(EFaceLocation::IMAX): {FInterpFaceI = iSurfIMAX; break;}
		case(EFaceLocation::JMIN): {FInterpFaceI = iSurfJMIN; break;}
		case(EFaceLocation::JMAX): {FInterpFaceI = iSurfJMAX; break;}
		default: ERROR("Face is unknown.");
	}

	// Return the function pointer.
	return FInterpFaceI;
}

//-----------------------------------------------------------------------------------

auto IInterface::GetFuncPointerInterpFaceJ
(
 void
)
 /*
	* Function which returns a function pointer for the interpolation on 
	* the (matching) face belonging to the j-marker.
	*/
{
	// Create a function pointer for the jface in the jmarker.
	std::function<void(const size_t, const as3double*, as3double*, as3double*, as3double*)> FInterpFaceJ;

	// Definitions of the four interpolation functions on the (matching) jface, 
	// which are given as lambda's that bind to std::function.
	auto jSurfIMIN = [this](auto... in){ mJTensorProductContainer->SurfaceIMIN(in...); };
	auto jSurfIMAX = [this](auto... in){ mJTensorProductContainer->SurfaceIMAX(in...); };
	auto jSurfJMIN = [this](auto... in){ mJTensorProductContainer->SurfaceJMIN(in...); };
	auto jSurfJMAX = [this](auto... in){ mJTensorProductContainer->SurfaceJMAX(in...); };

	// Assign the appropriate interpolation functions for the (matching) jface.
	switch(mJFace)
	{
		case(EFaceLocation::IMIN): {FInterpFaceJ = jSurfIMIN; break;}
		case(EFaceLocation::IMAX): {FInterpFaceJ = jSurfIMAX; break;}
		case(EFaceLocation::JMIN): {FInterpFaceJ = jSurfJMIN; break;}
		case(EFaceLocation::JMAX): {FInterpFaceJ = jSurfJMAX; break;}
		default: ERROR("Face is unknown.");
	}

	// Return the function pointer.
	return FInterpFaceJ;
}

//-----------------------------------------------------------------------------------

auto IInterface::GetFuncPointerResidualFaceI
(
 void
)
 /*
	* Function which returns a function pointer for the residual computation on 
	* the (owner) face belonging to the i-marker.
	*/
{
	// Create a function pointer for the iface in the imarker.
	std::function<void(const size_t, const as3double*, const as3double*, const as3double*, as3double*)> FResFaceI;

	// Definitions of the four interpolation functions on the (owner) iface, 
	// which are given as lambda's that bind to std::function.
	auto iSurfIMIN = [this](auto... in){ mITensorProductContainer->ResidualSurfaceIMIN(in...); };
	auto iSurfIMAX = [this](auto... in){ mITensorProductContainer->ResidualSurfaceIMAX(in...); };
	auto iSurfJMIN = [this](auto... in){ mITensorProductContainer->ResidualSurfaceJMIN(in...); };
	auto iSurfJMAX = [this](auto... in){ mITensorProductContainer->ResidualSurfaceJMAX(in...); };

	// Assign the appropriate interpolation functions for the (owner) iface.
	switch(mIFace)
	{
		case(EFaceLocation::IMIN): {FResFaceI = iSurfIMIN; break;}
		case(EFaceLocation::IMAX): {FResFaceI = iSurfIMAX; break;}
		case(EFaceLocation::JMIN): {FResFaceI = iSurfJMIN; break;}
		case(EFaceLocation::JMAX): {FResFaceI = iSurfJMAX; break;}
		default: ERROR("Face is unknown.");
	}

	// Return the function pointer.
	return FResFaceI;
}

//-----------------------------------------------------------------------------------

auto IInterface::GetFuncPointerResidualFaceJ
(
 void
)
 /*
	* Function which returns a function pointer for the residual computation on 
	* the (matching) face belonging to the j-marker.
	*/
{
	// Create a function pointer for the jface in the jmarker.
	std::function<void(const size_t, const as3double*, const as3double*, const as3double*, as3double*)> FResFaceJ;

	// Definitions of the four interpolation functions on the (matching) jface, 
	// which are given as lambda's that bind to std::function.
	auto jSurfIMIN = [this](auto... in){ mJTensorProductContainer->ResidualSurfaceIMIN(in...); };
	auto jSurfIMAX = [this](auto... in){ mJTensorProductContainer->ResidualSurfaceIMAX(in...); };
	auto jSurfJMIN = [this](auto... in){ mJTensorProductContainer->ResidualSurfaceJMIN(in...); };
	auto jSurfJMAX = [this](auto... in){ mJTensorProductContainer->ResidualSurfaceJMAX(in...); };

	// Assign the appropriate interpolation functions for the (matching) jface.
	switch(mJFace)
	{
		case(EFaceLocation::IMIN): {FResFaceJ = jSurfIMIN; break;}
		case(EFaceLocation::IMAX): {FResFaceJ = jSurfIMAX; break;}
		case(EFaceLocation::JMIN): {FResFaceJ = jSurfJMIN; break;}
		case(EFaceLocation::JMAX): {FResFaceJ = jSurfJMAX; break;}
		default: ERROR("Face is unknown.");
	}

	// Return the function pointer.
	return FResFaceJ;
}





