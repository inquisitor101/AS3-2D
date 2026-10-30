/*!
 * @file
 * @brief Regression checks for supported element dimensions and quadrature sizing.
 */

#include "compile_time_dimensions.hpp"
#include "integration_rule_structure.hpp"
#include "standard_element_structure.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>


// Ensure the metadata can be queried during constant evaluation.
static_assert(CSupportedElementDimensions::GetMaxSol1D() == 10);
static_assert(CSupportedElementDimensions::GetMaxInt1D() == 14);
static_assert(CSupportedElementDimensions::IsSupported(6, 7));
static_assert(CSupportedElementDimensions::IsSupported(6, 7, 4));
static_assert(!CSupportedElementDimensions::IsSupported(5, 6));
static_assert(!CSupportedElementDimensions::IsSupported(6, 7, 5));


namespace NTest_ElementDimensions
{
	void CheckQuadrature(ETypeQuadrature typeQuadrature,
	                     std::size_t     nPolyIntegrand,
	                     std::size_t     nExpected)
	 /*
		* Function that checks point selection against generated polynomial moments.
		*/
	{
		if( CIntegrationRule::ComputePointCount1D(typeQuadrature, nPolyIntegrand) != nExpected )
			throw std::runtime_error("Incorrect quadrature point count.");

		const CIntegrationRule rule(typeQuadrature, nPolyIntegrand);
		if( rule.GetnInt1D() != nExpected )
			throw std::runtime_error("Point selection and quadrature initialization disagree.");

		const auto &nodes   = rule.GetrInt1D();
		const auto &weights = rule.GetwInt1D();
		for(std::size_t degree=0; degree<=rule.GetnPolyExact1D(); degree++)
		{
			as3double integral = 0.0;
			for(std::size_t i=0; i<nExpected; i++)
			{
				integral += weights[i] * std::pow(nodes[i], static_cast<int>(degree));
			}
			const as3double exact = degree%2 ? 0.0 : 2.0/static_cast<as3double>(degree + 1);
			if( std::abs(integral - exact) > 1.0e-12 )
				throw std::runtime_error("Quadrature polynomial exactness changed.");
		}
	}
}


int main(int argc, char **argv)
 /*
	* Entry point for quadrature checks or an isolated standard-element construction.
	*/
{
	try
	{
		if( argc == 1 )
		{
			for(std::size_t degree=0; degree<=1; degree++)
			{
				NTest_ElementDimensions::CheckQuadrature(ETypeQuadrature::GAUSS_LEGENDRE, degree, 1);
				NTest_ElementDimensions::CheckQuadrature(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, degree, 2);
			}
			NTest_ElementDimensions::CheckQuadrature(ETypeQuadrature::GAUSS_LEGENDRE, 9, 5);
			NTest_ElementDimensions::CheckQuadrature(ETypeQuadrature::GAUSS_LEGENDRE, 10, 6);
			NTest_ElementDimensions::CheckQuadrature(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 9, 6);
			NTest_ElementDimensions::CheckQuadrature(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 10, 7);
			std::cout << "Element metadata and quadrature checks passed." << std::endl;
		}
		else
		{
			if( argc != 5 ) throw std::runtime_error("Expected: p bonus GL|LGL NODAL_DG|DG_SEM.");

			const auto nPoly  = static_cast<std::size_t>(std::stoull(argv[1]));
			const auto nBonus = static_cast<std::size_t>(std::stoull(argv[2]));
			const auto typeQuadrature = std::string(argv[3]) == "GL" ?
			                           ETypeQuadrature::GAUSS_LEGENDRE : ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO;
			const auto typeDGMethod = std::string(argv[4]) == "NODAL_DG" ?
			                         ETypeDGMethod::NODAL_DG : ETypeDGMethod::DG_SEM;
			const CStandardElement element(ETypeDOF::LGL, nPoly, nBonus, typeQuadrature, typeDGMethod);
			std::cout << "K=" << element.GetnSol1D() << " M=" << element.GetnInt1D() << std::endl;
		}

#ifdef ENABLE_NAN_CHECK
		NError::CheckFloatingError();
#endif
	}
	catch(const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}

	return 0;
}
