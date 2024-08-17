

//-----------------------------------------------------------------------------------
// Implementation of the templated/inlined functions in CPhysicalElement.
//-----------------------------------------------------------------------------------


void CPhysicalElement::ConvertGradParamToCartVolInt
(
 CWorkMatrixAS3<as3double> &dVarDr,
 CWorkMatrixAS3<as3double> &dVarDs
)
 /*
	* Function that converts the parametric gradient into a Cartesian 
	* gradient at the volume integration points. 
	*/
{
	// Determine inner and outer loops.
	const size_t nVar = dVarDr.row();
	const size_t nInt = dVarDr.col();

	// Convert the gradient from parametric to Cartesian coordinates.
	for(size_t i=0; i<nVar; i++)
	{
		for(size_t l=0; l<nInt; l++)
		{
			// Parametric derivatives.
			const as3double dudr = dVarDr(i,l);
			const as3double duds = dVarDs(i,l);
			
			// Metrics of transformation.
			const as3double drdx = mMetricInt2D(1,l);
			const as3double drdy = mMetricInt2D(2,l);
			const as3double dsdx = mMetricInt2D(3,l);
			const as3double dsdy = mMetricInt2D(4,l);

			// Compute the Cartesian derivatives w.r.t. x.
			dVarDr(i,l) = dudr*drdx + duds*dsdx;
			
			// Compute the Cartesian derivatives w.r.t. y.
			dVarDs(i,l) = dudr*drdy + duds*dsdy;
		}
	}
}


