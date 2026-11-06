#include "interface/IInterface.hpp"

#include <cstddef>

#include "compute_kernels/TensorProducts.hpp"
#include "geometry/face/CInterfaceFaceBatch.hpp"
#include "riemann_solver/IRiemannSolver.hpp"
#include "solver/CMultizoneSolver.hpp"
#include "solver/ISinglezoneSolver.hpp"
#include "element/CStandardElement.hpp"


//-----------------------------------------------------------------------------------
// IInterface member functions.
//-----------------------------------------------------------------------------------


IInterface::IInterface
(
 const CConfig               *config_container,
 const CMultizoneSolver      *multizone_solver_container,
 const CInterfaceFaceBatch   &interface_batch
)
 /*
  * Constructor that records paired-zone metadata and creates interface operators.
  */
{
  // Assign the general interface properties.
  mIName = interface_batch.GetOwnerFaceName();
  mJName = interface_batch.GetNeighborFaceName();

  mIZone = interface_batch.GetOwnerZoneIndex();
  mJZone = interface_batch.GetNeighborZoneIndex();

  mIFace = interface_batch.GetOwnerFaceLocation();
  mJFace = interface_batch.GetNeighborFaceLocation();

  mNFace = interface_batch.GetnFace(); 

	mNPolyI = multizone_solver_container->GetSinglezoneSolver(mIZone)->GetStandardElement()->GetnPolySol();
	mNPolyJ = multizone_solver_container->GetSinglezoneSolver(mJZone)->GetStandardElement()->GetnPolySol();
}

//-----------------------------------------------------------------------------------

IInterface::~IInterface
(
 void
)
 /*
	* Destructor for IInterface.
	*/
{

}
