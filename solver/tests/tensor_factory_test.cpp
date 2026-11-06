/*!
 * @file
 * @brief Regression checks for collocation detection and tensor factory selection.
 */

#include "CGenericFactory.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>


namespace NTest_TensorFactory
{
	static_assert(std::is_same_v<decltype(std::declval<CIntegrationRule>().GetTypeNodalDistribution()),
	                             decltype(std::declval<CStandardElement>().GetTypeDOFsSol())>);

	void Check(bool condition, const std::string &message)
	 /*
		* Function that reports a failed factory selection check.
		*/
	{
		if( !condition ) throw std::runtime_error(message);
	}

	template<std::size_t K>
	void CheckCollocated(void)
	 /*
		* Function that checks DG-SEM dispatch for one supported point count.
		*/
	{
		CStandardElement element(ETypeNodalDistribution::LGL, K-1, 0,
		                         ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, ETypeDGMethod::DG_SEM);
		Check(element.GetIsCollocated(), "DG-SEM nodes must be collocated.");
		auto tensor = CGenericFactory::CreateTensorContainer(&element, 4, ETypeDGMethod::DG_SEM);
		Check(dynamic_cast<CCollocatedTensorProduct<K,4>*>(tensor.get()) != nullptr,
		      "DG-SEM did not select the collocated specialization.");
	}
}


int main(int argc, char **argv)
 /*
	* Entry point for factory selection checks and invalid DG-SEM combinations.
	*/
{
	using namespace NTest_TensorFactory;
	try
	{
		// Quadrature and solution nodes use the same distribution type.
		const CIntegrationRule glRule(ETypeQuadrature::GAUSS_LEGENDRE, 9);
		const CIntegrationRule lglRule(ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, 9);
		Check(glRule.GetTypeNodalDistribution() == ETypeNodalDistribution::GL,
		      "Gauss-Legendre quadrature did not report GL nodes.");
		Check(lglRule.GetTypeNodalDistribution() == ETypeNodalDistribution::LGL,
		      "Gauss-Lobatto quadrature did not report LGL nodes.");
		Check(MapTypeNodalDistribution.at("GL") == ETypeNodalDistribution::GL &&
		      MapTypeNodalDistribution.at("LGL") == ETypeNodalDistribution::LGL &&
		      MapTypeNodalDistribution.at("EQD") == ETypeNodalDistribution::EQD,
		      "Incorrect common nodal distribution mapping.");

		// Representing GL nodes in the common enum does not enable GL solution face kernels.
		if( argc > 1 && std::string(argv[1]) == "gl_solution" )
		{
			CStandardElement element(ETypeNodalDistribution::GL, 5, 0,
			                         ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
			return 1;
		}

		// Check the complete configuration-to-factory selection without writing simulation output.
		if( argc == 3 && std::string(argv[1]) == "config" )
		{
			CConfig config(argv[2]);
			CStandardElement element(ETypeNodalDistribution::LGL, 5, config.GetBonusQuadratureFactor(),
			                         config.GetTypeQuadrature(), config.GetTypeDGMethod());
			auto tensor = CGenericFactory::CreateTensorContainer(&element, 4, config.GetTypeDGMethod());
			if( config.GetTypeDGMethod() == ETypeDGMethod::DG_SEM )
				Check(dynamic_cast<CCollocatedTensorProduct<6,4>*>(tensor.get()) != nullptr,
				      "Configured DG-SEM did not select collocated kernels.");
			else
				Check(dynamic_cast<CTensorProduct<6,6,4>*>(tensor.get()) != nullptr,
				      "Configured nodal DG did not select regular kernels.");
			std::cout << "Configured tensor selection passed." << std::endl;
			return 0;
		}

		// Equal counts with different nodal distributions must use the regular tensor.
		CStandardElement gl(ETypeNodalDistribution::LGL, 5, 0,
		                    ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
		Check(!gl.GetIsCollocated(), "LGL solution nodes are not collocated with GL quadrature.");
		if( argc > 1 && std::string(argv[1]) == "invalid" )
		{
			CGenericFactory::CreateTensorContainer(&gl, 4, ETypeDGMethod::DG_SEM);
			return 1;
		}
		auto regular = CGenericFactory::CreateTensorContainer(&gl, 4, ETypeDGMethod::NODAL_DG);
		Check(dynamic_cast<CTensorProduct<6,6,4>*>(regular.get()) != nullptr,
		      "Equal counts alone selected a collocated tensor.");

		// Matching LGL distributions with unequal counts are also noncollocated.
		CStandardElement lgl(ETypeNodalDistribution::LGL, 5, 0,
		                     ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, ETypeDGMethod::NODAL_DG);
		Check(!lgl.GetIsCollocated(), "Different LGL point counts must not be collocated.");
		regular = CGenericFactory::CreateTensorContainer(&lgl, 4, ETypeDGMethod::NODAL_DG);
		Check(dynamic_cast<CTensorProduct<6,7,4>*>(regular.get()) != nullptr,
		      "Noncollocated LGL nodal DG did not select the regular tensor.");

		// Equidistant solution nodes are not the LGL distribution.
		CStandardElement eqd(ETypeNodalDistribution::EQD, 5, 0,
		                     ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
		Check(!eqd.GetIsCollocated(), "Equidistant solution nodes must not be marked LGL collocated.");
		regular = CGenericFactory::CreateTensorContainer(&eqd, 4, ETypeDGMethod::NODAL_DG);
		Check(dynamic_cast<CTensorProduct<6,6,4>*>(regular.get()) != nullptr,
		      "Equidistant nodal DG did not select the regular tensor.");

		CheckCollocated<2>(); CheckCollocated<3>(); CheckCollocated<4>();
		CheckCollocated<5>(); CheckCollocated<6>(); CheckCollocated<7>();
		CheckCollocated<8>(); CheckCollocated<9>(); CheckCollocated<10>();
		std::cout << "All collocation and tensor factory checks passed." << std::endl;
	}
	catch(const std::exception &error)
	{
		std::cerr << error.what() << std::endl;
		return 1;
	}
}
