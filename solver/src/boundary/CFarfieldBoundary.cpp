#include "boundary/CFarfieldBoundary.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"

#include <cstddef>

#include "BoundaryParameters.hpp"


//-----------------------------------------------------------------------------------
// CFarfieldBoundary member functions.
//-----------------------------------------------------------------------------------


CFarfieldBoundary::CFarfieldBoundary
(
 const CConfig              *config_container,
 const CSinglezoneGeometry  *singlezone_geometry_container,
 const ISinglezoneSolver    *singlezone_solver_container,
 const CBoundaryFaceBatch   &boundary_batch,
 const IBoundaryParamMarker *param_bc
)
  : IBoundary(config_container,
              singlezone_geometry_container,
              singlezone_solver_container,
              boundary_batch,
              param_bc)
 /*
	* Constructor that initializes a farfield boundary.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }

  // Convert after validating the boundary condition type.
  const auto* farfield_param =
    static_cast<const CFarfieldParamMarker*>(param_bc);

  // Extract and assign the relevant prescribed values as the conservative state.
  mPrescribedDensity     = farfield_param->mDensity;
  mPrescribedMomentum[0] = mPrescribedDensity * farfield_param->mVelocity[0];
  mPrescribedMomentum[1] = mPrescribedDensity * farfield_param->mVelocity[1];
  mPrescribedTotalEnergy = C_OVGM1 *   farfield_param->mPressure
                         +  C_HALF * ( farfield_param->mVelocity[0] * mPrescribedMomentum[0] 
                                     + farfield_param->mVelocity[1] * mPrescribedMomentum[1] );

}

//-----------------------------------------------------------------------------------

CFarfieldBoundary::~CFarfieldBoundary
(
 void
)
 /*
	* Destructor for CFarfieldBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CFarfieldBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that imposes the farfield state from the prescribed values.
  */
{
  for(std::size_t l=0; l<met.col(); l++)
  {
    // Assemble the boundary state.
    var_b(0,l) = mPrescribedDensity;
    var_b(1,l) = mPrescribedMomentum[0];
    var_b(2,l) = mPrescribedMomentum[1];
    var_b(3,l) = mPrescribedTotalEnergy; 
  }
}
