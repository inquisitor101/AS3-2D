#include "solver/ISinglezoneSolver.hpp"

#include "compute_kernels/TensorProducts.hpp"
#include "CConfig.hpp"
#include "CGenericFactory.hpp"
#include "geometry/CSinglezoneGeometry.hpp"
#include "element/CPhysicalElement.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "element/CStandardElement.hpp"

#include <cstddef>


//-----------------------------------------------------------------------------------
// ISinglezoneSolver member functions.
//-----------------------------------------------------------------------------------


ISinglezoneSolver::ISinglezoneSolver
(
 const CConfig             *config_container,
 const CSinglezoneGeometry *singlezone_geometry_container,
 std::size_t                nVar
)
	:
		mZoneID( singlezone_geometry_container->GetZoneID() ),
		mTypeDGFormulation( config_container->GetTypeDGFormulation() )
 /*
	* Constructor that creates the zone reference element, tensor and Riemann solvers.
	*/
{
	// Instantiate a standard element object.
	mStandardElementContainer = CGenericFactory::CreateStandardElement(config_container,
                                                                     singlezone_geometry_container);

	// Instantiate a tensor-product object.
	mTensorProductContainer = CGenericFactory::CreateTensorContainer(mStandardElementContainer.get(),
                                                                   nVar, config_container->GetTypeDGMethod());

	// Instantiate a Riemann solver object.
	mRiemannSolverContainer = CGenericFactory::CreateRiemannSolverContainer(config_container,
			                                                                    config_container->GetTypeRiemannSolver(mZoneID));
}

//-----------------------------------------------------------------------------------

ISinglezoneSolver::~ISinglezoneSolver
(
 void
)
 /*
	* Destructor for ISinglezoneSolver.
	*/
{

}
