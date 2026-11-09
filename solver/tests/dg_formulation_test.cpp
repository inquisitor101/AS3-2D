/*!
 * @file
 * @brief Cartesian strong derivatives, conservation and free-stream checks on curved elements and zones.
 */

#include "CConfig.hpp"
#include "CGenericFactory.hpp"
#include "CSupportedElementDimensions.hpp"
#include "CIteration.hpp"
#include "element/CPhysicalElement.hpp"
#include "geometry/CElementGeometry.hpp"
#include "geometry/CMultizoneGeometry.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "parallel/COpenMP.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "solver/CEESinglezoneSolver.hpp"
#include "solver/CMultizoneSolver.hpp"
#include "temporal/CSSPRK3Temporal.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace NTest_DGFormulation
{
	void CheckClose(as3double actual, as3double expected, const char *message, as3double tolerance=3.e-10)
	 /*
		* Function that checks a scaled residual or solution difference.
		*/
	{
		if( !std::isfinite(actual) || !std::isfinite(expected) ||
		    std::abs(actual-expected) > tolerance*(C_ONE + std::abs(expected)) )
			throw std::runtime_error(std::string(message) + " Actual=" + std::to_string(actual) + " expected=" + std::to_string(expected));
	}

	void SetSolution(CPhysicalElement *element, const CMatrixAS3<as3double> &coord, bool constant)
	 /*
		* Function that assigns positive Euler states, including nonlinear and discontinuous traces.
		*/
	{
		auto& sol = element->mSol2D;
		for(std::size_t q=0; q<sol.col(); q++)
		{
			const as3double perturbation = constant ? C_ZERO : std::sin(coord(0,q) + 0.7*coord(1,q));
			const as3double rho = 1.0 + 0.1*perturbation;
			const as3double u   = 0.3 + 0.05*perturbation;
			const as3double v   = -0.2 + 0.03*perturbation;
			const as3double p   = 1.0 + 0.08*perturbation;
			sol(0,q) = rho;
			sol(1,q) = rho*u;
			sol(2,q) = rho*v;
			sol(3,q) = p/C_GM1 + C_HALF*rho*(u*u + v*v);
		}
	}


	CMatrixAS3<as3double> ComputeCartesianReference(ISinglezoneSolver *solver, CPhysicalElement *element, bool dense=false)
	 /*
		* Function that independently assembles the Cartesian strong residual using reference tensor derivatives.
		*/
	{
		const auto* standard = solver->GetStandardElement();
		const std::size_t K = standard->GetnInt1D(), Q = standard->GetnInt2D(), P = standard->GetnSol2D();
		CPoolMatrixAS3<as3double> work(2*4*P + 4*4*Q + 2*4*K);
		auto fx = work.GetWorkMatrixAS3(4, P);
		auto fy = work.GetWorkMatrixAS3(4, P);
		auto dxFx = work.GetWorkMatrixAS3(4, Q);
		auto dyFx = work.GetWorkMatrixAS3(4, Q);
		auto dxFy = work.GetWorkMatrixAS3(4, Q);
		auto dyFy = work.GetWorkMatrixAS3(4, Q);
		const auto& sol = element->mSol2D;
		for(std::size_t q=0; q<P; q++)
		{
			const as3double rho = sol(0,q), mx = sol(1,q), my = sol(2,q), energy = sol(3,q);
			const as3double pressure = C_GM1*(energy - C_HALF*(mx*mx + my*my)/rho);
			fx(0,q) = mx;
			fx(1,q) = mx*mx/rho + pressure;
			fx(2,q) = mx*my/rho;
			fx(3,q) = mx/rho*(energy + pressure);
			fy(0,q) = my;
			fy(1,q) = mx*my/rho;
			fy(2,q) = my*my/rho + pressure;
			fy(3,q) = my/rho*(energy + pressure);
		}
		auto* tensor = solver->GetTensorProduct();
		tensor->CompileTimeVolume(fx.data(), nullptr, dxFx.data(), dyFx.data());
		tensor->CompileTimeVolume(fy.data(), nullptr, dxFy.data(), dyFy.data());
		element->ConvertGradParamToCartVolInt(dxFx, dyFx);
		element->ConvertGradParamToCartVolInt(dxFy, dyFy);
		CMatrixAS3<as3double> weighted(4, Q), residual(4, P);
		for(std::size_t v=0; v<4; v++)
			for(std::size_t q=0; q<Q; q++)
				weighted(v,q) = -standard->GetwInt2D()[q]*element->mMetricInt2D(0,q)*(dxFx(v,q) + dyFy(v,q));

		if( dense )
		{
			for(std::size_t v=0; v<4; v++)
			{
				for(std::size_t q=0; q<Q; q++)
				{
					as3double divergence = C_ZERO;
					for(std::size_t l=0; l<P; l++)
					{
						divergence += element->mDerLagrangeInt2Dx(q,l)*fx(v,l) + element->mDerLagrangeInt2Dy(q,l)*fy(v,l);
					}
					weighted(v,q) = -standard->GetwInt2D()[q]*element->mMetricInt2D(0,q)*divergence;
				}
			}
		}

		tensor->CompileTimeResidualVolume(weighted.data(), nullptr, nullptr, residual.data());

		// Trace the same nodal flux polynomial used in the volume to each integration face.
		using Eval = void (ITensorProduct::*)(const as3double*, as3double*, as3double*, as3double*);
		using Lift = void (ITensorProduct::*)(const as3double*, const as3double*, const as3double*, as3double*);
		const Eval eval[] = {&ITensorProduct::CompileTimeSurfaceIMIN, &ITensorProduct::CompileTimeSurfaceIMAX,
		                     &ITensorProduct::CompileTimeSurfaceJMIN, &ITensorProduct::CompileTimeSurfaceJMAX};
		const Lift lift[] = {&ITensorProduct::CompileTimeResidualSurfaceIMIN, &ITensorProduct::CompileTimeResidualSurfaceIMAX,
		                     &ITensorProduct::CompileTimeResidualSurfaceJMIN, &ITensorProduct::CompileTimeResidualSurfaceJMAX};
		auto trace = work.GetWorkMatrixAS3(4, K);
		auto flux = work.GetWorkMatrixAS3(4, K);
		for(std::size_t face=0; face<4; face++)
		{
			(tensor->*eval[face])(fx.data(), trace.data(), nullptr, nullptr);
			(tensor->*eval[face])(fy.data(), flux.data(), nullptr, nullptr);
			const auto& met = element->GetSurfaceMetricInt(static_cast<EFaceLocation>(face));
			for(std::size_t v=0; v<4; v++)
			{
				for(std::size_t k=0; k<K; k++)
				{
					trace(v,k) = standard->GetwInt1D()[k]*met(0,k)*(met(1,k)*trace(v,k) + met(2,k)*flux(v,k));
				}
			}
			(tensor->*lift[face])(trace.data(), nullptr, nullptr, residual.data());
		}
		return residual;
	}



	void CheckInverseMass(CPhysicalElement &element, ITensorProduct *tensor, const CStandardElement *standard)
	 /*
		* Function that checks exact collocated diagonal storage and inversion of the quadrature mass action.
		*/
	{
		const std::size_t nSol2D = standard->GetnSol2D(), nInt2D = standard->GetnInt2D();
		const auto& weights = standard->GetwInt2D();
		if( standard->GetIsCollocated() )
		{
			if( element.mInvMassMatrix.row() != nSol2D || element.mInvMassMatrix.col() != 1 )
				throw std::runtime_error("Collocated inverse mass must use packed diagonal storage.");
			for(std::size_t q=0; q<nSol2D; q++)
			{
				const as3double expected = C_ONE/(element.mMetricInt2D(0,q)*weights[q]);
				if( element.mInvMassMatrix[q] != expected )
					throw std::runtime_error("Collocated inverse mass must contain exact reciprocals.");
			}
		}

		// Integrate each nodal field against the basis, then recover it by applying the inverse mass.
		// This also exercises non-collocated grids, including equal counts with different node distributions.
		CMatrixAS3<as3double> field(4, nSol2D), weighted(4, nInt2D), residual(4, nSol2D), scratch(1, nSol2D);
		for(std::size_t q=0; q<field.size(); q++) field[q] = 0.3 + std::sin(0.17*static_cast<as3double>(q));
		tensor->CompileTimeVolume(field.data(), weighted.data(), nullptr, nullptr);
		for(std::size_t v=0; v<4; v++)
		{
			for(std::size_t q=0; q<nInt2D; q++) weighted(v,q) *= weights[q]*element.mMetricInt2D(0,q);
		}
		tensor->CompileTimeResidualVolume(weighted.data(), nullptr, nullptr, residual.data());
		tensor->CompileTimeApplyInverseMassMatrix(element.mInvMassMatrix.data(), residual.data(),
		    standard->GetIsCollocated() ? nullptr : scratch.data());
		for(std::size_t q=0; q<field.size(); q++) CheckClose(residual[q], field[q], "Inverse mass did not recover the nodal field.");
	}

	void CheckNoncollocatedMatrices(CConfig &strong_config)
	 /*
		* Function that checks Cartesian basis derivatives with distinct solution and integration grids.
		*/
	{
		struct Case { ETypeQuadrature mQuadrature; std::size_t mBonus; };
		const Case cases[] = {{ETypeQuadrature::GAUSS_LEGENDRE, 0},
		                      {ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 0},
		                      {ETypeQuadrature::GAUSS_LEGENDRE, 1}};
		for(const auto& test : cases)
		{
			CStandardElement standard(ETypeNodalDistribution::LGL, 5, test.mBonus,
			                          test.mQuadrature, ETypeDGMethod::NODAL_DG);
			auto tensor = CGenericFactory::CreateTensorContainer(&standard, 4, ETypeDGMethod::NODAL_DG);
			CElementGeometry geometry(standard.GetnSol2D());
			auto& coord = geometry.GetCoordSolDOFs();
			const auto& rSol = standard.GetrSol1D();
			const auto& rInt = standard.GetrInt1D();
			const std::size_t K = standard.GetnSol1D(), M = standard.GetnInt1D();
			CMatrixAS3<as3double> polynomial(1, standard.GetnSol2D());
			for(std::size_t j=0; j<K; j++)
			{
				for(std::size_t i=0; i<K; i++)
				{
					const as3double r = rSol[i], s = rSol[j];
					coord(0,j*K+i) = r + 0.2*s + 0.08*r*s;
					coord(1,j*K+i) = 0.1*r + 0.9*s + 0.06*r*s;
					polynomial[j*K+i] = r*r*s + 0.3*r*s*s;
				}
			}

			// Build the element operators independently of the solver residual assembly.
			CPhysicalElement element(&strong_config, &standard, tensor.get(), &geometry, 0, 4);
			CheckInverseMass(element, tensor.get(), &standard);
			if( element.mDerLagrangeInt2Dx.row() != M*M || element.mDerLagrangeInt2Dx.col() != K*K ||
			    element.mDerLagrangeInt2Dy.row() != M*M || element.mDerLagrangeInt2Dy.col() != K*K )
				throw std::runtime_error("Incorrect non-collocated derivative matrix dimensions.");

			for(std::size_t j=0; j<M; j++)
			{
				for(std::size_t i=0; i<M; i++)
				{
					const std::size_t q = j*M+i;
					const as3double r = rInt[i], s = rInt[j];
					const as3double dPDr = 2.0*r*s + 0.3*s*s;
					const as3double dPDs = r*r + 0.6*r*s;
					const auto& met = element.mMetricInt2D;
					as3double dx = C_ZERO, dy = C_ZERO;
					for(std::size_t l=0; l<K*K; l++)
					{
						dx += element.mDerLagrangeInt2Dx(q,l)*polynomial[l];
						dy += element.mDerLagrangeInt2Dy(q,l)*polynomial[l];
					}
					CheckClose(dx, dPDr*met(1,q) + dPDs*met(3,q), "Non-collocated x-derivative mismatch.");
					CheckClose(dy, dPDr*met(2,q) + dPDs*met(4,q), "Non-collocated y-derivative mismatch.");
				}
			}
		}
		std::cout << "Non-collocated Cartesian and inverse-mass checks passed for GL, LGL and overintegration.\n";
	}

	void CheckLocal(CConfig &weak_config, CConfig &strong_config)
	 /*
		* Function that checks every supported collocated order on skew and curved mappings.
		*/
	{
		for(std::size_t p=1; p<=9; p++)
			for(bool curved : {false, true})
			{
				CSinglezoneGeometry zone("test", 0, p, 1, 1, ETypeNodalDistribution::LGL, !curved);
				CEESinglezoneSolver weak(&weak_config, &zone), strong(&strong_config, &zone);
				auto* geometry = zone.GetElementGeometry(0);
				auto& coord = geometry->GetCoordSolDOFs();
				const auto& nodes = weak.GetStandardElement()->GetrSol1D();
				const std::size_t K = p+1;
				for(std::size_t j=0; j<K; j++)
					for(std::size_t i=0; i<K; i++)
					{
						const as3double r = nodes[i], s = nodes[j];
						coord(0,j*K+i) = r + 0.2*s + (curved ? 0.08*r*s + 0.03*s*s : C_ZERO);
						coord(1,j*K+i) = 0.1*r + 0.9*s + (curved ? 0.06*r*s + 0.04*r*r : C_ZERO);
					}
				for(auto* solver : {&weak, &strong})
				{
					solver->GetPhysicalElement().push_back(std::make_unique<CPhysicalElement>(
					    solver == &weak ? &weak_config : &strong_config, const_cast<CStandardElement*>(solver->GetStandardElement()),
					    solver->GetTensorProduct(), geometry, 0, 4));
					CheckInverseMass(*solver->GetPhysicalElement(0), solver->GetTensorProduct(), solver->GetStandardElement());
				}

				const auto* element = strong.GetPhysicalElement(0);
				if( weak.GetPhysicalElement(0)->mDerLagrangeInt2Dx.size() ||
				    weak.GetPhysicalElement(0)->mDerLagrangeInt2Dy.size() )
					throw std::runtime_error("Weak elements must not allocate Cartesian derivative matrices.");
				for(std::size_t q=0; q<K*K; q++)
				{
					as3double dxX=C_ZERO, dyX=C_ZERO, dxY=C_ZERO, dyY=C_ZERO;
					for(std::size_t l=0; l<K*K; l++)
					{
						dxX += element->mDerLagrangeInt2Dx(q,l)*coord(0,l);
						dyX += element->mDerLagrangeInt2Dy(q,l)*coord(0,l);
						dxY += element->mDerLagrangeInt2Dx(q,l)*coord(1,l);
						dyY += element->mDerLagrangeInt2Dy(q,l)*coord(1,l);
					}
					CheckClose(dxX, C_ONE, "Incorrect dx/dx."); CheckClose(dyX, C_ZERO, "Incorrect dx/dy.");
					CheckClose(dxY, C_ZERO, "Incorrect dy/dx."); CheckClose(dyY, C_ONE, "Incorrect dy/dy.");
				}

				CPoolMatrixAS3<as3double> work(3*4*K*K);
				for(bool constant : {false, true})
				{
					for(auto* solver : {&weak, &strong}) SetSolution(solver->GetPhysicalElement(0), coord, constant);
					weak.ComputeElementLocalResidual(&zone, work, C_ZERO, 0);
					strong.ComputeElementLocalResidual(&zone, work, C_ZERO, 0);
					auto* element = strong.GetPhysicalElement(0);
					const auto& reference = weak.GetPhysicalElement(0)->mRes2D;
					const auto cartesian = ComputeCartesianReference(&strong, element);
					const auto dense = ComputeCartesianReference(&strong, element, true);
					for(std::size_t q=0; q<reference.size(); q++)
					{
						CheckClose(element->mRes2D[q], cartesian[q], "Cartesian strong local residual mismatch.");
						CheckClose(element->mRes2D[q], dense[q], "Factored/dense Cartesian residual mismatch.");
						// Variable curved metrics need not commute with differentiation of nonlinear fluxes.
						if( !curved || constant ) CheckClose(element->mRes2D[q], reference[q], "Affine/constant weak/strong mismatch.");
					}
					// The strong local volume and physical face terms have zero net integral.
					for(std::size_t v=0; v<4; v++)
					{
						as3double sum = C_ZERO;
						for(std::size_t q=0; q<K*K; q++) sum += element->mRes2D(v,q);
						CheckClose(sum, C_ZERO, "Local conservation failed.");
					}

					if( !constant ) continue;
					// Consistent Riemann fluxes must cancel all four local physical boundary fluxes.
					using Eval = void (ITensorProduct::*)(const as3double*, as3double*, as3double*, as3double*);
					using Lift = void (ITensorProduct::*)(const as3double*, const as3double*, const as3double*, as3double*);
					const Eval eval[] = {&ITensorProduct::CompileTimeSurfaceIMIN, &ITensorProduct::CompileTimeSurfaceIMAX,
					                     &ITensorProduct::CompileTimeSurfaceJMIN, &ITensorProduct::CompileTimeSurfaceJMAX};
					const Lift lift[] = {&ITensorProduct::CompileTimeResidualSurfaceIMIN, &ITensorProduct::CompileTimeResidualSurfaceIMAX,
					                     &ITensorProduct::CompileTimeResidualSurfaceJMIN, &ITensorProduct::CompileTimeResidualSurfaceJMAX};
					auto trace = work.GetWorkMatrixAS3(4, K);
					auto flux = work.GetWorkMatrixAS3(4, K);
					auto* tensor = strong.GetTensorProduct();
					for(std::size_t face=0; face<4; face++)
					{
						(tensor->*eval[face])(element->mSol2D.data(), trace.data(), nullptr, nullptr);
						strong.GetRiemannSolver()->ComputeFlux(strong.GetStandardElement()->GetwInt1D(),
						    element->GetSurfaceMetricInt(static_cast<EFaceLocation>(face)), trace, trace, flux);
						(tensor->*lift[face])(flux.data(), nullptr, nullptr, element->mRes2D.data());
					}
					for(std::size_t q=0; q<element->mRes2D.size(); q++)
						CheckClose(element->mRes2D[q], C_ZERO, "Curved-element free stream was not preserved.");
				}
			}
		std::cout << "Cartesian derivatives, inverse mass, local residual, conservation and free-stream checks passed for p=1..9.\n";
	}


	void CheckNoncollocatedResidual(CConfig &config, std::size_t p)
	 /*
		* Function that compares dense non-collocated residuals with independent tensor derivatives and lifting.
		*/
	{
		CSinglezoneGeometry zone("test", 0, p, 1, 1, ETypeNodalDistribution::LGL, false);
		CEESinglezoneSolver solver(&config, &zone);
		const auto* standard = solver.GetStandardElement();
		auto* geometry = zone.GetElementGeometry(0);
		auto& coord = geometry->GetCoordSolDOFs();
		const auto& nodes = standard->GetrSol1D();
		const std::size_t K = standard->GetnSol1D();
		for(std::size_t j=0; j<K; j++)
		{
			for(std::size_t i=0; i<K; i++)
			{
				const auto r=nodes[i], s=nodes[j];
				coord(0,j*K+i)=r+0.2*s+0.08*r*s+0.03*s*s;
				coord(1,j*K+i)=0.1*r+0.9*s+0.06*r*s+0.04*r*r;
			}
		}
		solver.GetPhysicalElement().push_back(std::make_unique<CPhysicalElement>(&config,
		    const_cast<CStandardElement*>(standard), solver.GetTensorProduct(), geometry, 0, 4));
		auto* element=solver.GetPhysicalElement(0);
		CPoolMatrixAS3<as3double> work(4*(2*standard->GetnSol2D()+standard->GetnInt2D()+2*standard->GetnInt1D()));
		for(bool constant : {false,true})
		{
			SetSolution(element,coord,constant);
			solver.ComputeElementLocalResidual(&zone,work,C_ZERO,0);
			const auto expected=ComputeCartesianReference(&solver,element);
			for(std::size_t q=0;q<expected.size();q++) CheckClose(element->mRes2D[q],expected[q],"Non-collocated strong residual mismatch.");
			for(std::size_t v=0;v<4;v++)
			{
				as3double sum=C_ZERO;
				for(std::size_t q=0;q<standard->GetnSol2D();q++) sum+=element->mRes2D(v,q);
				CheckClose(sum,C_ZERO,"Non-collocated local conservation failed.");
			}
		}
		std::cout << "Non-collocated local residual and conservation checks passed.\n";
	}

	void CheckMultizone(CConfig &weak_config, CConfig &strong_config)
	 /*
		* Function that checks face gathering, inverse mass, OpenMP traversal and all SSP-RK3 stages.
		*/
	{
		COpenMP openmp(&weak_config);
		CMultizoneGeometry geometry(&weak_config);
		CMultizoneSolver weak(&weak_config, &geometry), strong(&strong_config, &geometry);
		geometry.InitializeGridTopology(&openmp);
		for(auto* solvers : {&weak, &strong})
		{
			for(std::size_t z=0; z<geometry.GetnZone(); z++)
				solvers->GetSinglezoneSolver(z)->InitPhysicalElements(solvers == &weak ? &weak_config : &strong_config, &geometry);
			solvers->InitializeInterfaces(&weak_config, &geometry);
			solvers->InitializeBoundaries(&weak_config, &geometry);
		}
		CIteration weak_iteration(&weak_config, &weak), strong_iteration(&strong_config, &strong);
		for(bool constant : {false, true})
		{
			for(std::size_t z=0; z<geometry.GetnZone(); z++)
				for(std::size_t e=0; e<geometry.GetSinglezoneGeometry(z)->GetnElem(); e++)
					for(auto* solvers : {&weak, &strong})
					{
						auto* element = solvers->GetSinglezoneSolver(z)->GetPhysicalElement(e);
						SetSolution(element, geometry.GetSinglezoneGeometry(z)->GetElementGeometry(e)->GetCoordSolDOFs(), constant);
						if( !constant ) element->mSol2D(3,0) += 0.02;
					}
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
			{
				weak_iteration.GridSweep(&weak_config, &geometry, &openmp, &weak, C_ZERO);
				strong_iteration.GridSweep(&strong_config, &geometry, &openmp, &strong, C_ZERO);
			}
			as3double integral[4] = {}, scale[4] = {};
			as3double maxStrong = C_ZERO, maxWeak = C_ZERO;
			for(std::size_t z=0; z<geometry.GetnZone(); z++)
				for(std::size_t e=0; e<geometry.GetSinglezoneGeometry(z)->GetnElem(); e++)
				{
					auto* solver = strong.GetSinglezoneSolver(z);
					auto* actual = solver->GetPhysicalElement(e);
					const auto& expected = weak.GetSinglezoneSolver(z)->GetPhysicalElement(e)->mRes2D;
					const auto& weights = strong.GetSinglezoneSolver(z)->GetStandardElement()->GetwInt2D();
					auto cartesian = ComputeCartesianReference(solver, actual);
					for(const auto& face : actual->mResSurface)
						for(std::size_t q=0; q<cartesian.size(); q++) cartesian[q] += face[q];
					CMatrixAS3<as3double> scratch(1, weights.size());
					solver->GetTensorProduct()->CompileTimeApplyInverseMassMatrix(actual->mInvMassMatrix.data(), cartesian.data(), scratch.data());

					for(std::size_t v=0; v<4; v++)
						for(std::size_t q=0; q<actual->mSol2D.col(); q++)
						{
							CheckClose(actual->mRes2D(v,q), cartesian(v,q), "Multizone Cartesian strong mismatch.");
							if( constant )
							{
								// The imported fixture has a small constant-state residual, also present in weak form.
								// Bound the fixture error and compare the constant-state residual with weak form.
								CheckClose(actual->mRes2D(v,q), expected(v,q), "Constant-state weak/strong mismatch.");
								CheckClose(actual->mRes2D(v,q), C_ZERO, "Imported-grid free-stream tolerance exceeded.", 3.e-8);
								maxStrong = std::max(maxStrong, std::abs(actual->mRes2D(v,q)));
								maxWeak = std::max(maxWeak, std::abs(expected(v,q)));
							}
						}
					// Integrate the residual polynomial at quadrature points; these need not be solution nodes.
					CMatrixAS3<as3double> integrationResidual(4, weights.size());
					solver->GetTensorProduct()->CompileTimeVolume(actual->mRes2D.data(), integrationResidual.data(), nullptr, nullptr);
					for(std::size_t v=0; v<4; v++)
					{
						for(std::size_t q=0; q<weights.size(); q++)
						{
							const as3double value = weights[q]*actual->mMetricInt2D(0,q)*integrationResidual(v,q);
							integral[v] += value;
							scale[v] += std::abs(value);
						}
					}
				}
			for(std::size_t v=0; v<4; v++) CheckClose(integral[v]/(C_ONE+scale[v]), C_ZERO, "Periodic conservation failed.");

			if( constant ) std::cout << "Imported-grid free-stream maximum: strong=" << maxStrong << " weak=" << maxWeak << '\n';

			CSSPRK3Temporal temporal(&weak_config);
			for(std::size_t stage=0; stage<3; stage++)
			{
				weak_iteration.AdvanceStage(&geometry, &weak, &temporal, C_ZERO, 1.e-4, stage);
				strong_iteration.AdvanceStage(&geometry, &strong, &temporal, C_ZERO, 1.e-4, stage);
			}
			for(std::size_t z=0; z<geometry.GetnZone(); z++)
				for(std::size_t e=0; e<geometry.GetSinglezoneGeometry(z)->GetnElem(); e++)
				{
					const auto& actual = strong.GetSinglezoneSolver(z)->GetPhysicalElement(e)->mSol2D;
					const auto& expected = weak.GetSinglezoneSolver(z)->GetPhysicalElement(e)->mSol2D;
					for(std::size_t q=0; q<actual.size(); q++)
					{
						if( !std::isfinite(actual[q]) ) throw std::runtime_error("Nonfinite SSP-RK3 solution.");
						if( constant ) CheckClose(actual[q], expected[q], "Constant-state SSP-RK3 mismatch.");
					}
				}
		}
		std::cout << "Multizone Cartesian residual, conservation, free-stream and SSP-RK3 checks passed.\n";
	}
}

int main(int argc, char **argv)
 /*
	* Entry point for formulation regressions and configuration error paths.
	*/
{
	try
	{
		if( argc != 3 ) throw std::runtime_error("Expected weak and strong configuration filenames.");
		CConfig weak(argv[1]), strong(argv[2]);
		if( weak.GetTypeDGFormulation() != ETypeDGFormulation::WEAK_FORM ||
		    strong.GetTypeDGFormulation() != ETypeDGFormulation::STRONG_FORM )
			throw std::runtime_error("Incorrect DG formulation configuration.");
		NTest_DGFormulation::CheckNoncollocatedMatrices(strong);
		if( strong.GetTypeDGMethod() == ETypeDGMethod::DG_SEM ) NTest_DGFormulation::CheckLocal(weak, strong);
		else
		{
			for(std::size_t p=1; p<=9; p++)
			{
				CStandardElement standard(ETypeNodalDistribution::LGL, p, strong.GetBonusQuadratureFactor(),
				                          strong.GetTypeQuadrature(), ETypeDGMethod::NODAL_DG);
				if( CSupportedElementDimensions::IsSupported(standard.GetnSol1D(), standard.GetnInt1D(), 4) )
					NTest_DGFormulation::CheckNoncollocatedResidual(strong, p);
			}
		}
		NTest_DGFormulation::CheckMultizone(weak, strong);
	}
	catch(const std::exception &error)
	{
		std::cerr << error.what() << '\n';
		return 1;
	}
}
