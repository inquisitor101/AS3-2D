#include "test_unit_standard_element.hpp"


namespace
{
	TEST_F(CTest_SE, LocationDOFs)
	{
		for(auto& d: mElement)
		{
			CTest_SE::CheckBasisLocation1D(d.get());
			CTest_SE::CheckQuadrature(d.get());
		}
	}


}
