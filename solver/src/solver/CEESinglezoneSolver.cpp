#include "solver/CEESinglezoneSolver.hpp"
#include "data/CPoolMatrixAS3.hpp"

#include "compute_kernels/TensorProducts.hpp"
#include "CGenericFactory.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "element/CPhysicalElement.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "element/CStandardElement.hpp"

#include <cstddef>
#include <iostream>


//-----------------------------------------------------------------------------------
// CEESinglezoneSolver member functions.
//-----------------------------------------------------------------------------------


CEESinglezoneSolver::CEESinglezoneSolver
(
 const CConfig             *config_container,
 const CSinglezoneGeometry *singlezone_geometry_container
)
	:
		ISinglezoneSolver(config_container, singlezone_geometry_container, CEESinglezoneSolver::mNVar)
 /*
	* Constructor that initializes the Euler solver for one zone.
	*/
{

}

//-----------------------------------------------------------------------------------

CEESinglezoneSolver::~CEESinglezoneSolver
(
 void
)
 /*
	* Destructor for CEESinglezoneSolver.
	*/
{

}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::InitPhysicalElements
(
 const CConfig            *config_container,
 const CMultizoneGeometry *multizone_geometry_container
)
 /*
	* Function that constructs physical element states and metrics for the zone.
	*/
{	
	// Report output.
	std::cout << "    physical elements.... "; 

	// Get a reference to the current zone.
	auto* zone = multizone_geometry_container->GetSinglezoneGeometry(mZoneID);

	// Allocate memory for the physical elements in this solver.
	mPhysicalElementContainer.resize( zone->GetnElem() );
	
	// Instantiate physical element objects. The reason this needs a factory is because we might have 
	// a specialized solver (e.g. EE), but with a buffer zone, such as a PML. Hence, the elements
	// require additional variables. For now, avoid an interface class, as this object is the most 
	// performance-reliant data class in the code. Perhaps it is better to keep a single class, but 
	// instantiate it differently, depending on the type of solver and buffer layer specification.
	for(std::size_t i=0; i<mPhysicalElementContainer.size(); i++)
	{
		mPhysicalElementContainer[i] = CGenericFactory::CreatePhysicalElement(config_container,
				                                                                  mStandardElementContainer.get(),
																																					mTensorProductContainer.get(),
				                                                                  zone->GetElementGeometry(i),
																																					mZoneID, mNVar);
	}

	// Report output.
	std::cout << "Done." << std::endl;
}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::ComputeElementLocalResidual
(
 CSinglezoneGeometry       *singlezone_geometry_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElem
)
 /*
	* Function that overwrites one element residual with its local volume and physical face terms.
	*/
{
	if( mTypeDGFormulation == ETypeDGFormulation::WEAK_FORM )
	{
		// Extract the number of integration points in 2D.
		const std::size_t nInt2D = mStandardElementContainer->GetnInt2D();
		// Extract the quadrature integration weights on the standard element.
		auto&  wInt2D = mStandardElementContainer->GetwInt2D();

		// Borrow memory for the solution and its gradient.
		CWorkMatrixAS3<as3double> var    = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
		CWorkMatrixAS3<as3double> dVarDx = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
		CWorkMatrixAS3<as3double> dVarDy = workarray.GetWorkMatrixAS3(mNVar, nInt2D);

		// Extract the relevant physical element.
		auto& physical_element = mPhysicalElementContainer[iElem];

		// Extract the metrics at the volume solution points.
		auto& jac = physical_element->mMetricInt2D;
		// Reference to the current element solution.
		auto& sol = physical_element->mSol2D;
		// Reference to the current element residual.
		auto& res = physical_element->mRes2D;

		// Compute the solution and its (parametric) gradient at the volume integration points.
		mTensorProductContainer->CompileTimeVolume(sol.data(),
				                                       var.data(), dVarDx.data(), dVarDy.data());

		// Convert the gradient from parametric to Cartesian coordinates.
		physical_element->ConvertGradParamToCartVolInt(dVarDx, dVarDy);


		// Loop over the integration points and compute the weighted terms.
#pragma omp simd
	  for(std::size_t l=0; l<nInt2D; l++)
		{
			// Compute the weighting factor: Jacobian multiplied by the integration weight.
			// Note, the residual is defined on the RHS:
			//  dUDt = Res(vol) - Res(surf), thus the vol terms are +ve.
			const as3double weight = jac(0,l)*wInt2D[l];

			// Assemble the relevant (weighted) metrics in the x- and y-derivatives.
			const as3double wdrdx = weight * jac(1,l);
			const as3double wdrdy = weight * jac(2,l);
			const as3double wdsdx = weight * jac(3,l);
			const as3double wdsdy = weight * jac(4,l);

		  // Compute the primitive variables.
		  const as3double ovrho = C_ONE / var(0,l);
		  const as3double u     = ovrho * var(1,l);
		  const as3double v     = ovrho * var(2,l);
		  const as3double p     = C_GM1 * ( var(3,l) - C_HALF * (u * var(1,l) + v * var(2,l)) );

		  // Compute the inviscid flux in the x-direction.
		  const as3double fx0 =       var(1,l);
		  const as3double fx1 = u *   var(1,l) + p;
		  const as3double fx2 = v *   var(1,l);
		  const as3double fx3 = u * ( var(3,l) + p );

		  // Compute the inviscid flux in the y-direction.
		  const as3double fy0 =       var(2,l);
		  const as3double fy1 = u *   var(2,l);
		  const as3double fy2 = v *   var(2,l) + p;
		  const as3double fy3 = v * ( var(3,l) + p );

			// Compute the residual flux terms in the r-parametric direction.
			dVarDx(0,l) = fx0 * wdrdx + fy0 * wdrdy;
			dVarDx(1,l) = fx1 * wdrdx + fy1 * wdrdy;
			dVarDx(2,l) = fx2 * wdrdx + fy2 * wdrdy;
			dVarDx(3,l) = fx3 * wdrdx + fy3 * wdrdy;

			// Compute the residual flux terms in the s-parametric direction.
			dVarDy(0,l) = fx0 * wdsdx + fy0 * wdsdy;
			dVarDy(1,l) = fx1 * wdsdx + fy1 * wdsdy;
			dVarDy(2,l) = fx2 * wdsdx + fy2 * wdsdy;
			dVarDy(3,l) = fx3 * wdsdx + fy3 * wdsdy;
		}

		// Scatter the volume residual terms back to the actual residual. Note, this resets the residual!
		mTensorProductContainer->CompileTimeResidualVolume(nullptr,
	                                                     dVarDx.data(),
	                                                     dVarDy.data(),
	                                                     res.data());
	}
	else if( mTypeDGFormulation == ETypeDGFormulation::STRONG_FORM &&
	         mStandardElementContainer->GetIsCollocated() )
	{
		// Extract the collocated integration-point counts and quadrature weights.
		const std::size_t nInt1D = mStandardElementContainer->GetnInt1D();
		const std::size_t nInt2D = mStandardElementContainer->GetnInt2D();
		const auto& wInt1D       = mStandardElementContainer->GetwInt1D();
		const auto& wInt2D       = mStandardElementContainer->GetwInt2D();

		// Extract the solution, metrics and precomputed Cartesian differentiation matrices.
		auto& physical_element = mPhysicalElementContainer[iElem];
		const auto& sol        = physical_element->mSol2D;
		const auto& jac        = physical_element->mMetricInt2D;
		auto& res              = physical_element->mRes2D;

		// Borrow memory for the Cartesian fluxes. The collocated solution needs no interpolation.
		CWorkMatrixAS3<as3double> fluxX = workarray.GetWorkMatrixAS3(mNVar, nInt2D);
		CWorkMatrixAS3<as3double> fluxY = workarray.GetWorkMatrixAS3(mNVar, nInt2D);

		// Compute the Cartesian fluxes at every integration point before assembling their derivatives.
#pragma omp simd
		for(std::size_t l=0; l<nInt2D; l++)
		{
			// Compute the primitive variables.
			const as3double ovrho = C_ONE / sol(0,l);
			const as3double u     = ovrho * sol(1,l);
			const as3double v     = ovrho * sol(2,l);
			const as3double p     = C_GM1 * ( sol(3,l) - C_HALF * (u * sol(1,l) + v * sol(2,l)) );

			// Compute the inviscid flux in the x-direction.
			fluxX(0,l) =      sol(1,l);
			fluxX(1,l) = u *  sol(1,l) + p;
			fluxX(2,l) = v *  sol(1,l);
			fluxX(3,l) = u * (sol(3,l) + p);

			// Compute the inviscid flux in the y-direction.
			fluxY(0,l) =      sol(2,l);
			fluxY(1,l) = u *  sol(2,l);
			fluxY(2,l) = v *  sol(2,l) + p;
			fluxY(3,l) = v * (sol(3,l) + p);
		}


		// Extract the shared 1D reference differentiation matrix: dL1D(i,k) = L_k'(r_i).
		// The same matrix applies in the s-direction, since both directions use the same reference nodes.
		const auto& dL1D = mStandardElementContainer->GetDerLagrangeInt1D();

		// Compute the divergence of the Cartesian fluxes for each solution variable.
		for(std::size_t l=0; l<mNVar; l++)
		{
			// Loop over the evaluation point (r_i, s_j). Collocation makes this a solution point as well.
			for(std::size_t j=0; j<nInt1D; j++)
			{
				for(std::size_t i=0; i<nInt1D; i++)
				{
					// Flatten the evaluation-point index, with the i-index varying fastest.
					const std::size_t q = j*nInt1D + i;

					// Initialize the reference derivatives of both Cartesian flux components at q.
					as3double dFxDr = C_ZERO; as3double dFxDs = C_ZERO;
					as3double dFyDr = C_ZERO; as3double dFyDs = C_ZERO;

					// Under collocation, the transverse interpolation factor is the identity:
					//   dF/dr(i,j) = sum_k dL1D(i,k) * F(k,j), along row j;
					//   dF/ds(i,j) = sum_k dL1D(j,k) * F(i,k), along col i.
					// Thus, each derivative needs only nInt1D terms, rather than nInt2D terms.
#pragma omp simd reduction(+:dFxDr,dFxDs,dFyDr,dFyDs)
					for(std::size_t k=0; k<nInt1D; k++)
					{
						// Select the flux points (r_k, s_j) and (r_i, s_k), respectively.
						const std::size_t iRow = j*nInt1D + k;
						const std::size_t iCol = k*nInt1D + i;

						// Evaluate the derivatives of basis function k at r_i and s_j.
						const as3double dLDr = dL1D(i,k);
						const as3double dLDs = dL1D(j,k);

						// Accumulate the r-derivatives of Fx and Fy from the current row.
						dFxDr += dLDr*fluxX(l, iRow);
						dFyDr += dLDr*fluxY(l, iRow);

						// Accumulate the s-derivatives of Fx and Fy from the current column.
						dFxDs += dLDs*fluxX(l, iCol);
						dFyDs += dLDs*fluxY(l, iCol);
					}

					// Extract the inverse mapping metrics at the evaluation point q.
					const as3double drdx = jac(1,q);
					const as3double drdy = jac(2,q);
					const as3double dsdx = jac(3,q);
					const as3double dsdy = jac(4,q);

					// Apply the Cartesian chain rule, to compute the divergence of the fluxes.
					const as3double divergence = drdx * dFxDr + dsdx * dFxDs + drdy * dFyDr + dsdy * dFyDs;

					// Initialize the volume residual with -J*w*(dFx/dx + dFy/dy).
					// The negative sign follows from placing the residual on the RHS of the evolution equation.
					res(l,q) = -jac(0,q) * wInt2D[q] * divergence;
				}
			}
		}

		// Collocation places the face integration points directly on the volume grid.
		// Reuse the Cartesian fluxes at those points; no face interpolation is needed.
		const auto& metIMin = physical_element->mMetricIntIMin1D;
		const auto& metIMax = physical_element->mMetricIntIMax1D;
		const auto& metJMin = physical_element->mMetricIntJMin1D;
		const auto& metJMax = physical_element->mMetricIntJMax1D;

		// Each face contributes w*Js*(nx*Fx + ny*Fy), where Js is the surface Jacobian.
		// The stored unit normals already point outward, including on IMIN and JMIN.
		for(std::size_t v=0; v<mNVar; v++)
		{
			// IMIN: i = 0; integrate along the face by varying j.
#pragma omp simd
			for(std::size_t j=0; j<nInt1D; j++)
			{
				// Corresponding volume index, with i varying fastest in q = j*nInt1D + i.
				const std::size_t q = j*nInt1D;

				// Extract the outward unit normal and the quadrature-weighted surface Jacobian.
				const as3double nx     = metIMin(1,j);
				const as3double ny     = metIMin(2,j);
				const as3double weight = wInt1D[j]*metIMin(0,j);

				// Add the physical normal flux of this solution variable to its volume residual.
				const as3double normalFlux = nx*fluxX(v,q) + ny*fluxY(v,q);
				res(v,q) += weight*normalFlux;
			}

			// IMAX: i = nInt1D-1; integrate along the face by varying j.
#pragma omp simd
			for(std::size_t j=0; j<nInt1D; j++)
			{
				// Corresponding volume index, with i varying fastest in q = j*nInt1D + i.
				const std::size_t q = j*nInt1D + nInt1D-1;

				// Extract the outward unit normal and the quadrature-weighted surface Jacobian.
				const as3double nx     = metIMax(1,j);
				const as3double ny     = metIMax(2,j);
				const as3double weight = wInt1D[j]*metIMax(0,j);

				// Add the physical normal flux of this solution variable to its volume residual.
				const as3double normalFlux = nx*fluxX(v,q) + ny*fluxY(v,q);
				res(v,q) += weight*normalFlux;
			}

			// JMIN: j = 0; integrate along the face by varying i.
#pragma omp simd
			for(std::size_t i=0; i<nInt1D; i++)
			{
				// Corresponding volume index, with i varying fastest in q = j*nInt1D + i.
				const std::size_t q = i;

				// Extract the outward unit normal and the quadrature-weighted surface Jacobian.
				const as3double nx     = metJMin(1,i);
				const as3double ny     = metJMin(2,i);
				const as3double weight = wInt1D[i]*metJMin(0,i);

				// Add the physical normal flux of this solution variable to its volume residual.
				const as3double normalFlux = nx*fluxX(v,q) + ny*fluxY(v,q);
				res(v,q) += weight*normalFlux;
			}

			// JMAX: j = nInt1D-1; integrate along the face by varying i.
#pragma omp simd
			for(std::size_t i=0; i<nInt1D; i++)
			{
				// Corresponding volume index, with i varying fastest in q = j*nInt1D + i.
				const std::size_t q = (nInt1D-1)*nInt1D + i;

				// Extract the outward unit normal and the quadrature-weighted surface Jacobian.
				const as3double nx     = metJMax(1,i);
				const as3double ny     = metJMax(2,i);
				const as3double weight = wInt1D[i]*metJMax(0,i);

				// Add the physical normal flux of this solution variable to its volume residual.
				const as3double normalFlux = nx*fluxX(v,q) + ny*fluxY(v,q);
				res(v,q) += weight*normalFlux;
			}

		} // Corner points receive one contribution from each of their two adjoining faces.

	}
	else if( mTypeDGFormulation == ETypeDGFormulation::STRONG_FORM )
	{
		// Extract the number of solution points in 2D.
		const std::size_t nSol2D = mStandardElementContainer->GetnSol2D();
		// Extract the number of integration points in 1D and 2D.
		const std::size_t nInt1D = mStandardElementContainer->GetnInt1D();
		const std::size_t nInt2D = mStandardElementContainer->GetnInt2D();
		// Extract the quadrature integration weights on the standard element.
		const auto& wInt1D = mStandardElementContainer->GetwInt1D();
		const auto& wInt2D = mStandardElementContainer->GetwInt2D();

		// Extract the relevant physical element.
		auto& physical_element = mPhysicalElementContainer[iElem];
		// Reference to the current element solution and residual.
		const auto& sol = physical_element->mSol2D;
		auto&       res = physical_element->mRes2D;
		// Extract the metrics at the volume integration points.
		const auto& jac = physical_element->mMetricInt2D;
		// Extract the Cartesian basis derivatives at the integration points.
		// Rows correspond to integration points; columns correspond to solution basis functions.
		const auto& dLDx = physical_element->mDerLagrangeInt2Dx;
		const auto& dLDy = physical_element->mDerLagrangeInt2Dy;

		// Borrow memory for the Cartesian fluxes at the solution points.
		// These nodal values define the flux polynomials differentiated in the volume and traced to the faces.
		CWorkMatrixAS3<as3double> fluxX = workarray.GetWorkMatrixAS3(mNVar, nSol2D);
		CWorkMatrixAS3<as3double> fluxY = workarray.GetWorkMatrixAS3(mNVar, nSol2D);

		// Loop over the solution points and compute the Cartesian fluxes.
#pragma omp simd
		for(std::size_t l=0; l<nSol2D; l++)
		{
			// Compute the primitive variables.
			const as3double ovrho = C_ONE / sol(0,l);
			const as3double u     = ovrho * sol(1,l);
			const as3double v     = ovrho * sol(2,l);
			const as3double p     = C_GM1 * ( sol(3,l) - C_HALF * (u * sol(1,l) + v * sol(2,l)) );

			// Compute the inviscid flux in the x-direction.
			fluxX(0,l) =      sol(1,l);
			fluxX(1,l) = u *  sol(1,l) + p;
			fluxX(2,l) = v *  sol(1,l);
			fluxX(3,l) = u * (sol(3,l) + p);

			// Compute the inviscid flux in the y-direction.
			fluxY(0,l) =      sol(2,l);
			fluxY(1,l) = u *  sol(2,l);
			fluxY(2,l) = v *  sol(2,l) + p;
			fluxY(3,l) = v * (sol(3,l) + p);
		}

		// Keep the volume temporary in its own scope, so the same workspace can be borrowed for the faces later.
		{
			// Borrow memory for the weighted flux divergence at the volume integration points.
			CWorkMatrixAS3<as3double> weightedDivergence = workarray.GetWorkMatrixAS3(mNVar, nInt2D);

			// Compute the divergence of the Cartesian fluxes for each solution variable.
			for(std::size_t v=0; v<mNVar; v++)
			{
				// Loop over the integration points where the derivatives are evaluated.
				for(std::size_t q=0; q<nInt2D; q++)
				{
					as3double divergence = C_ZERO;

					// Sum over the solution-flux coefficients. The Cartesian metrics at q are already in dLDx/dLDy.
#pragma omp simd reduction(+:divergence)
					for(std::size_t l=0; l<nSol2D; l++)
					{
						divergence += dLDx(q,l) * fluxX(v,l) + dLDy(q,l) * fluxY(v,l);
					}

					// Weight the divergence by the quadrature weight and mapping Jacobian at q.
					// The residual is on the RHS, hence the negative volume-divergence sign.
					const as3double weight = -jac(0,q) * wInt2D[q];
					weightedDivergence(v,q) = weight * divergence;
				}
			}

			// Integrate the weighted divergence against the test basis and scatter it to the solution points.
			// The derivatives are already included above, so only value lifting is needed here.
			// Note, this resets the element residual!
			mTensorProductContainer->CompileTimeResidualVolume(weightedDivergence.data(),
			                                                   nullptr, nullptr, res.data());
		} // Release the volume temporary before borrowing the face buffers, respecting the pool's LIFO order.

		// Borrow memory for the Cartesian fluxes at the face integration points.
		// Reuse these buffers for all four faces, retaining the nodal volume fluxes for interpolation.
		CWorkMatrixAS3<as3double> faceFluxX = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
		CWorkMatrixAS3<as3double> faceFluxY = workarray.GetWorkMatrixAS3(mNVar, nInt1D);

		// Extract the surface metrics at the integration points of each face.
		const auto& metIMin = physical_element->mMetricIntIMin1D;
		const auto& metIMax = physical_element->mMetricIntIMax1D;
		const auto& metJMin = physical_element->mMetricIntJMin1D;
		const auto& metJMax = physical_element->mMetricIntJMax1D;

		// Interpolate the same flux polynomials used in the volume, rather than reevaluating F at an interpolated state.
		// For nonlinear fluxes, these operations differ; using the same polynomial keeps the local terms consistent.

		// IMIN: interpolate the Cartesian fluxes to this face's integration points.
		mTensorProductContainer->CompileTimeSurfaceIMIN(fluxX.data(), faceFluxX.data(), nullptr, nullptr);
		mTensorProductContainer->CompileTimeSurfaceIMIN(fluxY.data(), faceFluxY.data(), nullptr, nullptr);

		// Compute the weighted outward normal flux. The stored unit normal already includes its outward sign.
		for(std::size_t v=0; v<mNVar; v++)
		{
#pragma omp simd
			for(std::size_t k=0; k<nInt1D; k++)
			{
				const as3double nx     = metIMin(1,k);
				const as3double ny     = metIMin(2,k);
				const as3double weight = wInt1D[k] * metIMin(0,k);
				const as3double normalFlux = nx * faceFluxX(v,k) + ny * faceFluxY(v,k);

				// Reuse faceFluxX for the weighted normal flux after both Cartesian components have been read.
				faceFluxX(v,k) = weight * normalFlux;
			}
		}

		// Lift the IMIN physical face contribution into the element residual. This increments the residual.
		mTensorProductContainer->CompileTimeResidualSurfaceIMIN(faceFluxX.data(), nullptr, nullptr, res.data());

		// IMAX: interpolate the Cartesian fluxes to this face's integration points.
		mTensorProductContainer->CompileTimeSurfaceIMAX(fluxX.data(), faceFluxX.data(), nullptr, nullptr);
		mTensorProductContainer->CompileTimeSurfaceIMAX(fluxY.data(), faceFluxY.data(), nullptr, nullptr);

		// Compute the weighted outward normal flux. The stored unit normal already includes its outward sign.
		for(std::size_t v=0; v<mNVar; v++)
		{
#pragma omp simd
			for(std::size_t k=0; k<nInt1D; k++)
			{
				const as3double nx     = metIMax(1,k);
				const as3double ny     = metIMax(2,k);
				const as3double weight = wInt1D[k] * metIMax(0,k);
				const as3double normalFlux = nx * faceFluxX(v,k) + ny * faceFluxY(v,k);

				// Reuse faceFluxX for the weighted normal flux after both Cartesian components have been read.
				faceFluxX(v,k) = weight * normalFlux;
			}
		}

		// Lift the IMAX physical face contribution into the element residual. This increments the residual.
		mTensorProductContainer->CompileTimeResidualSurfaceIMAX(faceFluxX.data(), nullptr, nullptr, res.data());

		// JMIN: interpolate the Cartesian fluxes to this face's integration points.
		mTensorProductContainer->CompileTimeSurfaceJMIN(fluxX.data(), faceFluxX.data(), nullptr, nullptr);
		mTensorProductContainer->CompileTimeSurfaceJMIN(fluxY.data(), faceFluxY.data(), nullptr, nullptr);

		// Compute the weighted outward normal flux. The stored unit normal already includes its outward sign.
		for(std::size_t v=0; v<mNVar; v++)
		{
#pragma omp simd
			for(std::size_t k=0; k<nInt1D; k++)
			{
				const as3double nx     = metJMin(1,k);
				const as3double ny     = metJMin(2,k);
				const as3double weight = wInt1D[k] * metJMin(0,k);
				const as3double normalFlux = nx * faceFluxX(v,k) + ny * faceFluxY(v,k);

				// Reuse faceFluxX for the weighted normal flux after both Cartesian components have been read.
				faceFluxX(v,k) = weight * normalFlux;
			}
		}

		// Lift the JMIN physical face contribution into the element residual. This increments the residual.
		mTensorProductContainer->CompileTimeResidualSurfaceJMIN(faceFluxX.data(), nullptr, nullptr, res.data());

		// JMAX: interpolate the Cartesian fluxes to this face's integration points.
		mTensorProductContainer->CompileTimeSurfaceJMAX(fluxX.data(), faceFluxX.data(), nullptr, nullptr);
		mTensorProductContainer->CompileTimeSurfaceJMAX(fluxY.data(), faceFluxY.data(), nullptr, nullptr);

		// Compute the weighted outward normal flux. The stored unit normal already includes its outward sign.
		for(std::size_t v=0; v<mNVar; v++)
		{
#pragma omp simd
			for(std::size_t k=0; k<nInt1D; k++)
			{
				const as3double nx     = metJMax(1,k);
				const as3double ny     = metJMax(2,k);
				const as3double weight = wInt1D[k] * metJMax(0,k);
				const as3double normalFlux = nx * faceFluxX(v,k) + ny * faceFluxY(v,k);

				// Reuse faceFluxX for the weighted normal flux after both Cartesian components have been read.
				faceFluxX(v,k) = weight * normalFlux;
			}
		}

		// Lift the JMAX physical face contribution into the element residual. This increments the residual.
		mTensorProductContainer->CompileTimeResidualSurfaceJMAX(faceFluxX.data(), nullptr, nullptr, res.data());

	}
	else ERROR("Unknown DG formulation in element-local residual.");
}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::ComputeElementCoupledResidualIDir
(
 const CSinglezoneGeometry *singlezone_geometry_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElemL
)
 /*
	* Function that computes one same-zone I-face flux and fills both local face buffers.
	*/
{
  // Deduce the right element's index.
  const std::size_t iElemR = iElemL + 1;

	// Extract the number of integration points in 1D.
	std::size_t nInt1D = mStandardElementContainer->GetnInt1D();
	// Extract the quadrature integration weights in 1D on the standard element.
	auto&  wInt1D = mStandardElementContainer->GetwInt1D();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varL = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> varR = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, nInt1D);

	// Get a pointer to the respective left and right element, w.r.t. this face.
	auto& elemL = mPhysicalElementContainer[iElemL];
	auto& elemR = mPhysicalElementContainer[iElemR];

	// Extract the metrics at the integration points on the imax face of the left element.
	auto& metL = elemL->mMetricIntIMax1D;
	// Reference to the left element solution.
	auto& solL = elemL->mSol2D;
	// Reference to the right element residual (temporary).
	auto& resL = elemL->GetSurfaceResidual(EFaceLocation::IMAX);

	// Reset the values in the resL, before proceeding.
	for(std::size_t l=0; l<resL.size(); l++) resL[l] = C_ZERO;

	// Reference to the right element solution.
	auto& solR = elemR->mSol2D;
	// Reference to the right element residual.
	auto& resR = elemR->GetSurfaceResidual(EFaceLocation::IMIN);
  resR.reset();

	// Compute the solution on the integration nodes of the left element.
	mTensorProductContainer->CompileTimeSurfaceIMAX(solL.data(),
                                                  varL.data(), nullptr, nullptr); 

	// Compute the solution on the integration nodes of the right element.
	mTensorProductContainer->CompileTimeSurfaceIMIN(solR.data(),
		                                              varR.data(), nullptr, nullptr); 

	// Compute the flux state, weighted by the integration nodes and metrics. 
	// Notice, this is based on the left state, which is the outward-pointing
	// normal vector.
	mRiemannSolverContainer->ComputeFlux(wInt1D, metL, varL, varR, flux);

	// Compute the residual on the left  element, which is on the IMAX boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceIMAX(flux.data(), 
                                                          nullptr, nullptr, resL.data());

	// For local conservation, negate the flux, since it leaves the left element to enter the right element.
	for(std::size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

	// Compute the residual on the right element, which is on the IMIN boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceIMIN(flux.data(), 
                                                          nullptr, nullptr, resR.data());
}

//-----------------------------------------------------------------------------------

void CEESinglezoneSolver::ComputeElementCoupledResidualJDir
(
 const CSinglezoneGeometry *singlezone_geometry_container,
 CPoolMatrixAS3<as3double> &workarray,
 as3double                  localtime,
 std::size_t                iElemB
)
 /*
	* Function that computes one same-zone J-face flux and fills both local face buffers.
	*/
{
  // Extract the number of elements in the ith-direction in this zone.
	const std::size_t niElem = singlezone_geometry_container->GetniElem();	
  // Deduce the top element's index.
  const std::size_t iElemT = iElemB + niElem;

	// Extract the number of integration points in 1D.
	std::size_t nInt1D = mStandardElementContainer->GetnInt1D();
	// Extract the quadrature integration weights in 1D on the standard element.
	auto&  wInt1D = mStandardElementContainer->GetwInt1D();

	// Borrow memory for the solution on the two sides.
	CWorkMatrixAS3<as3double> varB = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> varT = workarray.GetWorkMatrixAS3(mNVar, nInt1D);
	CWorkMatrixAS3<as3double> flux = workarray.GetWorkMatrixAS3(mNVar, nInt1D);

	// Get a pointer to the respective left and right element, w.r.t. this face.
	auto& elemB = mPhysicalElementContainer[iElemB];
	auto& elemT = mPhysicalElementContainer[iElemT];

	// Extract the metrics at the integration points on the imax face of the left element.
	auto& metB = elemB->mMetricIntJMax1D;
	// Reference to the left element solution.
	auto& solB = elemB->mSol2D;
	// Reference to the right element residual (temporary).
	auto& resB = elemB->GetSurfaceResidual(EFaceLocation::JMAX);

	// Reset the values in the resL, before proceeding.
	for(std::size_t l=0; l<resB.size(); l++) resB[l] = C_ZERO;

	// Reference to the right element solution.
	auto& solT = elemT->mSol2D;
	// Reference to the right element residual.
	auto& resT = elemT->GetSurfaceResidual(EFaceLocation::JMIN);
  resT.reset();

	// Compute the solution on the integration nodes of the left element.
	mTensorProductContainer->CompileTimeSurfaceJMAX(solB.data(),
		                                              varB.data(), nullptr, nullptr); 

	// Compute the solution on the integration nodes of the right element.
	mTensorProductContainer->CompileTimeSurfaceJMIN(solT.data(),
		                                              varT.data(), nullptr, nullptr); 

	// Compute the flux state, weighted by the integration nodes and metrics. 
	// Notice, this is based on the left state, which is the outward-pointing
	// normal vector.
	mRiemannSolverContainer->ComputeFlux(wInt1D, metB, varB, varT, flux);

	// Compute the residual on the left  element, which is on the IMAX boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceJMAX(flux.data(), 
                                                          nullptr, nullptr, resB.data());

	// For local conservation, negate the flux, since it leaves the left element to enter the right element.
	for(std::size_t l=0; l<flux.size(); l++) flux[l] *= -C_ONE;

	// Compute the residual on the right element, which is on the IMIN boundary.
	mTensorProductContainer->CompileTimeResidualSurfaceJMIN(flux.data(), 
                                                          nullptr, nullptr, resT.data());
}
