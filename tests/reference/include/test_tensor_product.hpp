#pragma once

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "factory_structure.hpp"
#include "standard_element_structure.hpp"
#include "tensor_structure.hpp"
#include "test_linear_algebra.hpp"
#include "gtest/gtest.h"


class CTest_TP : public ::testing::Test
{
	public:
		void SetUp(void) override;
		void TearDown(void) override { } 

		std::unique_ptr<CConfig>                       mConfig  = nullptr;
		std::vector<std::unique_ptr<CStandardElement>> mElement;

		void CheckVolumeSol(CStandardElement *element);
		void CheckVolumeDerSolR(CStandardElement *element);
		void CheckVolumeDerSolS(CStandardElement *element);

		void CheckSurfaceSolIMIN(CStandardElement *element);
		void CheckSurfaceSolIMAX(CStandardElement *element);
		void CheckSurfaceSolJMIN(CStandardElement *element);
		void CheckSurfaceSolJMAX(CStandardElement *element);
		
		void CheckSurface_dSolDrIMIN(CStandardElement *element);
		void CheckSurface_dSolDrIMAX(CStandardElement *element);
		void CheckSurface_dSolDrJMIN(CStandardElement *element);
		void CheckSurface_dSolDrJMAX(CStandardElement *element);

		void CheckSurface_dSolDsIMIN(CStandardElement *element);
		void CheckSurface_dSolDsIMAX(CStandardElement *element);
		void CheckSurface_dSolDsJMIN(CStandardElement *element);
		void CheckSurface_dSolDsJMAX(CStandardElement *element);

		void CheckResidualVolumeSource(CStandardElement  *element);
		void CheckResidualVolumeDerSolR(CStandardElement *element);
		void CheckResidualVolumeDerSolS(CStandardElement *element);
		void CheckResidualVolumeTotal(CStandardElement   *element);


		static void KroneckerProduct(CMatrixAS3<as3double> &A, 
				                         CMatrixAS3<as3double> &B,
													       CMatrixAS3<as3double> &C);


};
