/*!
 * @file
 * @brief Inline implementations for element solution storage, metric transformations and mass matrices.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the templated/inlined functions in CPhysicalElement.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that computes the Cartesian gradient from the parametric gradient at the volume integration
 *   points.
 *
 * @param[in,out] dVarDr r-derivative on input and x-derivative on output.
 * @param[in,out] dVarDs s-derivative on input and y-derivative on output.
 */
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
	const std::size_t nVar = dVarDr.row();
	const std::size_t nInt = dVarDr.col();

	// Convert the gradient from parametric to Cartesian coordinates.
	for(std::size_t i=0; i<nVar; i++)
	{
		for(std::size_t l=0; l<nInt; l++)
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

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that computes the Cartesian gradient from the parametric gradient at the IMIN surface
 *   integration points.
 *
 * @param[in,out] dVarDr r-derivative on input and x-derivative on output.
 * @param[in,out] dVarDs s-derivative on input and y-derivative on output.
 */
void CPhysicalElement::ConvertGradParamToCartSurfIMinInt
(
 CWorkMatrixAS3<as3double> &dVarDr,
 CWorkMatrixAS3<as3double> &dVarDs
)
 /*
	* Function that converts the parametric gradient into a Cartesian 
	* gradient at the IMIN surface integration points. 
	*/
{
	// Determine inner and outer loops.
	const std::size_t nVar = dVarDr.row();
	const std::size_t nInt = dVarDr.col();

	// Convert the gradient from parametric to Cartesian coordinates.
	for(std::size_t i=0; i<nVar; i++)
	{
		for(std::size_t l=0; l<nInt; l++)
		{
			// Parametric derivatives.
			const as3double dudr = dVarDr(i,l);
			const as3double duds = dVarDs(i,l);
			
			// Metrics of transformation.
			const as3double drdx = mMetricIntIMin1D(1,l);
			const as3double drdy = mMetricIntIMin1D(2,l);
			const as3double dsdx = mMetricIntIMin1D(3,l);
			const as3double dsdy = mMetricIntIMin1D(4,l);

			// Compute the Cartesian derivatives w.r.t. x.
			dVarDr(i,l) = dudr*drdx + duds*dsdx;
			
			// Compute the Cartesian derivatives w.r.t. y.
			dVarDs(i,l) = dudr*drdy + duds*dsdy;
		}
	}
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that computes the Cartesian gradient from the parametric gradient at the IMAX surface
 *   integration points.
 *
 * @param[in,out] dVarDr r-derivative on input and x-derivative on output.
 * @param[in,out] dVarDs s-derivative on input and y-derivative on output.
 */
void CPhysicalElement::ConvertGradParamToCartSurfIMaxInt
(
 CWorkMatrixAS3<as3double> &dVarDr,
 CWorkMatrixAS3<as3double> &dVarDs
)
 /*
	* Function that converts the parametric gradient into a Cartesian 
	* gradient at the IMAX surface integration points. 
	*/
{
	// Determine inner and outer loops.
	const std::size_t nVar = dVarDr.row();
	const std::size_t nInt = dVarDr.col();

	// Convert the gradient from parametric to Cartesian coordinates.
	for(std::size_t i=0; i<nVar; i++)
	{
		for(std::size_t l=0; l<nInt; l++)
		{
			// Parametric derivatives.
			const as3double dudr = dVarDr(i,l);
			const as3double duds = dVarDs(i,l);
			
			// Metrics of transformation.
			const as3double drdx = mMetricIntIMax1D(1,l);
			const as3double drdy = mMetricIntIMax1D(2,l);
			const as3double dsdx = mMetricIntIMax1D(3,l);
			const as3double dsdy = mMetricIntIMax1D(4,l);

			// Compute the Cartesian derivatives w.r.t. x.
			dVarDr(i,l) = dudr*drdx + duds*dsdx;
			
			// Compute the Cartesian derivatives w.r.t. y.
			dVarDs(i,l) = dudr*drdy + duds*dsdy;
		}
	}
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that computes the Cartesian gradient from the parametric gradient at the JMIN surface
 *   integration points.
 *
 * @param[in,out] dVarDr r-derivative on input and x-derivative on output.
 * @param[in,out] dVarDs s-derivative on input and y-derivative on output.
 */
void CPhysicalElement::ConvertGradParamToCartSurfJMinInt
(
 CWorkMatrixAS3<as3double> &dVarDr,
 CWorkMatrixAS3<as3double> &dVarDs
)
 /*
	* Function that converts the parametric gradient into a Cartesian 
	* gradient at the JMIN surface integration points. 
	*/
{
	// Determine inner and outer loops.
	const std::size_t nVar = dVarDr.row();
	const std::size_t nInt = dVarDr.col();

	// Convert the gradient from parametric to Cartesian coordinates.
	for(std::size_t i=0; i<nVar; i++)
	{
		for(std::size_t l=0; l<nInt; l++)
		{
			// Parametric derivatives.
			const as3double dudr = dVarDr(i,l);
			const as3double duds = dVarDs(i,l);
			
			// Metrics of transformation.
			const as3double drdx = mMetricIntJMin1D(1,l);
			const as3double drdy = mMetricIntJMin1D(2,l);
			const as3double dsdx = mMetricIntJMin1D(3,l);
			const as3double dsdy = mMetricIntJMin1D(4,l);

			// Compute the Cartesian derivatives w.r.t. x.
			dVarDr(i,l) = dudr*drdx + duds*dsdx;
			
			// Compute the Cartesian derivatives w.r.t. y.
			dVarDs(i,l) = dudr*drdy + duds*dsdy;
		}
	}
}

//-----------------------------------------------------------------------------------

/*!
 * @brief Function that computes the Cartesian gradient from the parametric gradient at the JMAX surface
 *   integration points.
 *
 * @param[in,out] dVarDr r-derivative on input and x-derivative on output.
 * @param[in,out] dVarDs s-derivative on input and y-derivative on output.
 */
void CPhysicalElement::ConvertGradParamToCartSurfJMaxInt
(
 CWorkMatrixAS3<as3double> &dVarDr,
 CWorkMatrixAS3<as3double> &dVarDs
)
 /*
	* Function that converts the parametric gradient into a Cartesian 
	* gradient at the JMAX surface integration points. 
	*/
{
	// Determine inner and outer loops.
	const std::size_t nVar = dVarDr.row();
	const std::size_t nInt = dVarDr.col();

	// Convert the gradient from parametric to Cartesian coordinates.
	for(std::size_t i=0; i<nVar; i++)
	{
		for(std::size_t l=0; l<nInt; l++)
		{
			// Parametric derivatives.
			const as3double dudr = dVarDr(i,l);
			const as3double duds = dVarDs(i,l);
			
			// Metrics of transformation.
			const as3double drdx = mMetricIntJMax1D(1,l);
			const as3double drdy = mMetricIntJMax1D(2,l);
			const as3double dsdx = mMetricIntJMax1D(3,l);
			const as3double dsdy = mMetricIntJMax1D(4,l);

			// Compute the Cartesian derivatives w.r.t. x.
			dVarDr(i,l) = dudr*drdx + duds*dsdx;
			
			// Compute the Cartesian derivatives w.r.t. y.
			dVarDs(i,l) = dudr*drdy + duds*dsdy;
		}
	}
}


