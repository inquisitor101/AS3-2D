#include "boundary/CRiemannPressureOutletBoundary.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"

#include <cmath>
#include <cstddef>

#include "BoundaryParameters.hpp"


//-----------------------------------------------------------------------------------
// CRiemannPressureOutletBoundary member functions.
//-----------------------------------------------------------------------------------


CRiemannPressureOutletBoundary::CRiemannPressureOutletBoundary
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
	* Constructor that stores the prescribed outlet pressure.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }

  // Convert after validating the boundary condition type.
  const auto* outlet_param =
    static_cast<const CPressureOutletRiemannParamMarker*>(param_bc);

  // Extract and assign the relevant prescribed values.
  mPrescribedPressure = outlet_param->mPressure;
}

//-----------------------------------------------------------------------------------

CRiemannPressureOutletBoundary::~CRiemannPressureOutletBoundary
(
 void
)
 /*
	* Destructor for CRiemannPressureOutletBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CRiemannPressureOutletBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that forms the exterior outlet state with prescribed static pressure.
  */
{
  // Abbreviations.
  const as3double pb  = mPrescribedPressure;

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

    // Compute the Riemann invariant which is extrapolated.
    const as3double riemann = C_TWO * a * C_OVGM1 + un; 

    // Compute the specific entropy.
    const as3double s = p * std::pow( ovrho, C_GMA );

    // Compute the density, based on the extrapolated state.
    const as3double rb = std::pow( pb/s, C_OVGMA );

    // Compute the speed of sound, based on the extrapolated state.
    const as3double ab = std::sqrt( C_GMA * pb / rb );

    // Compute the normal velocity, based on the extrapolated state.
    const as3double ubn = riemann - C_TWO * ab * C_OVGM1; 
  
    // Compute the velocity of the extrapolated state.
    const as3double ub = u + (ubn - un) * nx;
    const as3double vb = v + (ubn - un) * ny;

    // Compute the specific kinetic energy of the extrapolated state.
    const as3double ekb = C_HALF * ( ub*ub + vb*vb );

    // Assemble the boundary state.
    var_b(0,l) = rb;
    var_b(1,l) = rb * ub;
    var_b(2,l) = rb * vb;
    var_b(3,l) = pb * C_OVGM1 + rb * ekb;
  }
}
