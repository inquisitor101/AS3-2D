#include "test_unit_tensor_product.hpp"


namespace
{

	TEST_F(CTest_TP, Volume)
	{
		for(auto& d: mElement)
		{
			CTest_TP::CheckVolumeSol(d.get());
			CTest_TP::CheckVolumeDerSolR(d.get());
			CTest_TP::CheckVolumeDerSolS(d.get());
		}
	}


	TEST_F(CTest_TP, Surface)
	{
		for(auto& d: mElement)
		{
			CTest_TP::CheckSurfaceSolIMIN(d.get());
			CTest_TP::CheckSurfaceSolIMAX(d.get());
			CTest_TP::CheckSurfaceSolJMIN(d.get());
			CTest_TP::CheckSurfaceSolJMAX(d.get());
			
			CTest_TP::CheckSurface_dSolDrIMIN(d.get());
			CTest_TP::CheckSurface_dSolDrIMAX(d.get());
			CTest_TP::CheckSurface_dSolDrJMIN(d.get());
			CTest_TP::CheckSurface_dSolDrJMAX(d.get());

			CTest_TP::CheckSurface_dSolDsIMIN(d.get());
			CTest_TP::CheckSurface_dSolDsIMAX(d.get());
			CTest_TP::CheckSurface_dSolDsJMIN(d.get());
			CTest_TP::CheckSurface_dSolDsJMAX(d.get());
		}
	}

	TEST_F(CTest_TP, ResidualVolume)
	{
		for(auto& d: mElement)
		{
			CTest_TP::CheckResidualVolumeSource(d.get());
			CTest_TP::CheckResidualVolumeDerSolR(d.get());
			CTest_TP::CheckResidualVolumeDerSolS(d.get());
			CTest_TP::CheckResidualVolumeTotal(d.get());
		}

	}

}
