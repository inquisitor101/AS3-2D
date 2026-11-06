#include "boundary/CRiemannStaticMassInletBoundary.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"

#include <cmath>
#include <cstddef>

#include "BoundaryParameters.hpp"


//-----------------------------------------------------------------------------------
// CRiemannStaticMassInletBoundary member functions.
//-----------------------------------------------------------------------------------


CRiemannStaticMassInletBoundary::CRiemannStaticMassInletBoundary
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
	* Constructor that stores the prescribed inlet density and velocity.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }

  // Convert after validating the boundary condition type.
  const auto* inlet_param =
    static_cast<const CStaticMassFlowInletRiemannParamMarker*>(param_bc);

  // Extract and assign the relevant prescribed values.
  mPrescribedDensity  = inlet_param->mDensity;
  mPrescribedVelocity = inlet_param->mVelocity;
}

//-----------------------------------------------------------------------------------

CRiemannStaticMassInletBoundary::~CRiemannStaticMassInletBoundary
(
 void
)
 /*
	* Destructor for CRiemannStaticMassInletBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CRiemannStaticMassInletBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that forms the exterior inlet state from prescribed density and velocity.
  */
{
  // Abbreviations.
  const as3double rb  = mPrescribedDensity;
  const as3double ub  = mPrescribedVelocity[0];
  const as3double vb  = mPrescribedVelocity[1];
  const as3double ekb = C_HALF * ( ub*ub + vb*vb );
  
  // Abbreviations involving gamma.
  const as3double gm1ovt  = C_HALF * C_GM1; 
  const as3double tgovgm1 = C_TWO * C_GMA * C_OVGM1;

  for(std::size_t l=0; l<met.col(); l++)
  {
    const as3double nx = met(1,l);
		const as3double ny = met(2,l);

    const as3double ovrho = C_ONE / var_e(0,l);
    const as3double u     = ovrho * var_e(1,l);
    const as3double v     = ovrho * var_e(2,l);
    const as3double p     = C_GM1  * (   var_e(3,l) 
                          - C_HALF * ( u*var_e(1,l) + v*var_e(2,l) ) );

    // Compute the acoustic information.
    const as3double a2 = C_GMA * p * ovrho;
    const as3double a  = std::sqrt( a2 );

    // Compute the normal velocity component.
    const as3double un = u * nx + v * ny;

    // Compute the normal velocity, based on the prescribed conditions.
    const as3double ubn = ub * nx + vb * ny; 
  
    // Compute the intermediate sound-speed ratio a_star/a along the interior
    // isentropic acoustic wave curve, using the prescribed normal velocity.
    const as3double ratio = C_ONE + gm1ovt * (un - ubn) / a;

    // Recover the corresponding static pressure. The boundary density is
    // prescribed independently when assembling the conservative state below.
    const as3double pb = p * std::pow(ratio, tgovgm1);

    // Assemble the boundary state.
    var_b(0,l) = rb;
    var_b(1,l) = rb * ub;
    var_b(2,l) = rb * vb;
    var_b(3,l) = pb * C_OVGM1 + rb * ekb;
  }
}
