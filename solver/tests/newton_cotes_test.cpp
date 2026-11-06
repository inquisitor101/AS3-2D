/*!
 * @file
 * @brief Regression checks for equidistant Newton-Cotes rules and their limits.
 */

#include "quadrature/CIntegrationRule.hpp"
#include "quadrature/legendre_gauss_lobatto.hpp"
#include "quadrature/newton_cotes.hpp"
#include "element/CStandardElement.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>


static_assert(NNewtonCotes::Validate());
static_assert(!NNewtonCotes::IsAvailable(0));
static_assert(NNewtonCotes::IsAvailable(8));
static_assert(!NNewtonCotes::IsAvailable(9));
static_assert(!NNewtonCotes::IsAvailable(10));


namespace NTest_NewtonCotes
{
	void Check(bool condition, const std::string &message)
	 /*
		* Function that reports a failed regression check.
		*/
	{
		if( !condition ) throw std::runtime_error(message);
	}

	void CheckMoments(std::span<const as3double> points,
	                  std::span<const as3double> weights,
	                  std::size_t                degreeExact)
	 /*
		* Function that checks the advertised exactness and the first inexact moment.
		*/
	{
		for(std::size_t degree=0; degree<=degreeExact+1; degree++)
		{
			as3double integral = 0.0;
			for(std::size_t i=0; i<points.size(); i++)
				integral += weights[i] * std::pow(points[i], static_cast<int>(degree));
			const as3double exact = degree%2 ? 0.0 : 2.0/static_cast<as3double>(degree + 1);
			if( degree <= degreeExact )
				Check(std::abs(integral - exact) < 1.0e-12, "Incorrect Newton-Cotes polynomial moment.");
			else
				Check(std::abs(integral - exact) > 1.0e-6, "Newton-Cotes exactness is understated.");
		}
	}

	void CheckTables(void)
	 /*
		* Function that exercises every stored rule, including even point counts.
		*/
	{
		for(std::size_t n=1; n<=8; n++)
		{
			const auto points  = NNewtonCotes::GetPoints(n);
			const auto weights = NNewtonCotes::GetWeights(n);
			Check(points.size() == n && weights.size() == n, "Incorrect table dimensions.");
			for(std::size_t i=0; i<n; i++)
			{
				Check(weights[i] > 0.0, "Nonpositive Newton-Cotes weight.");
				const as3double expected = n == 1 ? 0.0 : -1.0 + 2.0*static_cast<as3double>(i)
				                                                              /static_cast<as3double>(n-1);
				Check(std::abs(points[i] - expected) < 1.0e-14, "Nodes are not equidistant.");
			}
			CheckMoments(points, weights, n - 1 + n%2);
		}

		// Independent familiar rules: trapezoidal, Simpson and Simpson 3/8.
		Check(NNewtonCotes::GetWeights(2)[0] == 1.0, "Incorrect trapezoidal weight.");
		Check(NNewtonCotes::GetWeights(3)[0] == 1.0/3.0, "Incorrect Simpson endpoint weight.");
		Check(NNewtonCotes::GetWeights(3)[1] == 4.0/3.0, "Incorrect Simpson midpoint weight.");
		Check(NNewtonCotes::GetWeights(4)[0] == 0.25, "Incorrect Simpson 3/8 endpoint weight.");
		Check(NNewtonCotes::GetWeights(4)[1] == 0.75, "Incorrect Simpson 3/8 interior weight.");

		// A positive symmetric LGL rule has endpoints and moments, but wrong spacing.
		Check(!NIntegrationRuleValidation::ValidateRule(ETypeQuadrature::NEWTON_COTES,
		        NLegendreGaussLobatto::GetPoints(5), NLegendreGaussLobatto::GetWeights(5)),
		      "Accepted non-equidistant Newton-Cotes points.");
		const std::array<as3double, 3> points = {-1.0, 0.0, 1.0};
		const std::array<as3double, 3> weights = {1.5, -1.0, 1.5};
		Check(!NIntegrationRuleValidation::ValidateRule(ETypeQuadrature::NEWTON_COTES, points, weights),
		      "Accepted negative Newton-Cotes weights.");
	}

	void CheckIntegrationRules(void)
	 /*
		* Function that checks degree selection, metadata and tensor-product moments.
		*/
	{
		const auto type = ETypeQuadrature::NEWTON_COTES;
		Check(MapTypeQuadrature.at("NEWTON_COTES") == type, "Missing quadrature configuration keyword.");
		const std::array<std::size_t, 8> counts = {1, 1, 3, 3, 5, 5, 7, 7};
		for(std::size_t degree=0; degree<counts.size(); degree++)
		{
			const std::size_t n = counts[degree];
			const CIntegrationRule rule(type, degree);
			Check(CIntegrationRule::ComputePointCount1D(type, degree) == n, "Incorrect degree-to-count selection.");
			Check(rule.GetnInt1D() == n && rule.GetnInt2D() == n*n, "Incorrect generated dimensions.");
			Check(rule.GetTypeQuadrature() == type, "Incorrect quadrature metadata.");
			Check(rule.GetTypeNodalDistribution() == ETypeNodalDistribution::EQD, "Incorrect nodal distribution.");
			Check(rule.GetnPolyExact1D() == n, "Incorrect exactness metadata.");
			Check(rule.IsExactForPolynomial1D(degree), "Requested exactness was not provided.");
			Check(!rule.IsExactForPolynomial1D(n+1), "Overstated polynomial exactness.");

			const auto &points  = rule.GetrInt1D();
			const auto &weights = rule.GetwInt1D();
			CheckMoments({points.data(), n}, {weights.data(), n}, n);
			for(std::size_t a=0; a<=n; a++)
			{
				for(std::size_t b=0; b<=n; b++)
				{
					as3double integral = 0.0;
					for(std::size_t j=0; j<n; j++)
					{
						for(std::size_t i=0; i<n; i++)
						{
							const as3double w = rule.GetwInt2D()[j*n+i];
							Check(w > 0.0 && w == weights[i]*weights[j], "Incorrect tensor-product weight.");
							integral += w * std::pow(points[i], static_cast<int>(a))
							              * std::pow(points[j], static_cast<int>(b));
						}
					}
					const as3double exact = a%2 || b%2 ? 0.0 : 4.0/(static_cast<as3double>(a+1)
					                                                            *static_cast<as3double>(b+1));
					Check(std::abs(integral - exact) < 1.0e-12, "Incorrect tensor-product moment.");
				}
			}
			const CIntegrationRule sameRule(type, n);
			Check(rule == sameRule, "Equivalent Newton-Cotes rules compare unequal.");
		}

		// Exercise reference-element construction with equidistant integration nodes.
		const CStandardElement element(ETypeNodalDistribution::EQD, 3, 0, type, ETypeDGMethod::NODAL_DG);
		Check(element.GetnSol1D() == 4 && element.GetnInt1D() == 7, "Incorrect standard-element dimensions.");
	}
}


int main(int argc, char **argv)
 /*
	* Entry point for positive checks or isolated ERROR-path checks.
	*/
{
	try
	{
		if( argc == 3 )
		{
			const std::string mode = argv[1];
			const std::size_t value = std::string(argv[2]) == "max" ? std::numeric_limits<std::size_t>::max()
			                                                        : static_cast<std::size_t>(std::stoull(argv[2]));
			if( mode == "degree" )
			{
				const CIntegrationRule rule(ETypeQuadrature::NEWTON_COTES, value);
			}
			else if( mode == "count" )
				(void)CIntegrationRule::ComputePointCount1D(ETypeQuadrature::NEWTON_COTES, value);
			else if( mode == "points" ) (void)NNewtonCotes::GetPoints(value);
			else if( mode == "weights" ) (void)NNewtonCotes::GetWeights(value);
			else if( mode == "element" )
			{
				const CStandardElement element(ETypeNodalDistribution::EQD, value, 0,
				                               ETypeQuadrature::NEWTON_COTES, ETypeDGMethod::NODAL_DG);
			}
			else throw std::runtime_error("Unknown test mode.");
		}
		else
		{
			if( argc != 1 ) throw std::runtime_error("Expected: [degree|count|points|weights|element value].");
			NTest_NewtonCotes::CheckTables();
			NTest_NewtonCotes::CheckIntegrationRules();
			std::cout << "Newton-Cotes tables, exactness and tensor-product checks passed." << std::endl;
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
