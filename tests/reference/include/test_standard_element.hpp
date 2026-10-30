#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "standard_element_structure.hpp"
#include "test_linear_algebra.hpp"
#include "test_tensor_product.hpp"
#include "gtest/gtest.h"


class CTest_SE : public ::testing::Test
{
	public:
		void SetUp(void) override;
		void TearDown(void) override { } 

		std::unique_ptr<CConfig>                       mConfig  = nullptr;
		std::vector<std::unique_ptr<CStandardElement>> mElement;

		void CheckBasisLocation1D(CStandardElement *element);
		void CheckQuadrature(CStandardElement *element);


		as3double ComputeLagrangePolynomial(as3double               x,
				                                std::vector<as3double> &r,
																				size_t                  i);

		void ComputeLocationDOFs1D(size_t                 nDOFs1D, 
				                       ETypeDOF               typeDOFs, 
												       CMatrixAS3<as3double> &rDOFs1D);

		void ComputeQuadrature1D(size_t                   nInt1D,
				                     CMatrixAS3<as3double>   &rInt1D,
														 CMatrixAS3<as3double>   &wInt1D);
};
