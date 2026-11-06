#include "boundary/CSymmetryBoundary.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"

#include <cstddef>

#include "BoundaryParameters.hpp"


//-----------------------------------------------------------------------------------
// CSymmetryBoundary member functions.
//-----------------------------------------------------------------------------------


CSymmetryBoundary::CSymmetryBoundary
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
	* Constructor that initializes an inviscid symmetry boundary.
	*/
{
  // Consistency check.
  if( param_bc->GetTypeBC() != this->GetTypeBC() )
  {
    ERROR("Input boundary parameter has a wrong boundary condition.");
  }
}

//-----------------------------------------------------------------------------------

CSymmetryBoundary::~CSymmetryBoundary
(
 void
)
 /*
	* Destructor for CSymmetryBoundary.
	*/
{

}

//-----------------------------------------------------------------------------------

void CSymmetryBoundary::ComputeBoundaryState
(
 const CMatrixAS3<as3double>     &met,
 const CWorkMatrixAS3<as3double> &var_e,
 CWorkMatrixAS3<as3double>       &var_b
) const
 /*
  * Function that reflects normal velocity to construct the exterior symmetry state.
  */
{
  for(std::size_t l=0; l<met.col(); l++)
  {
    const as3double nx = met(1,l);
		const as3double ny = met(2,l);

    const as3double ovrho = C_ONE / var_e(0,l);
          as3double u     = ovrho * var_e(1,l);
          as3double v     = ovrho * var_e(2,l);
    const as3double p     = C_GM1  * (   var_e(3,l) 
                          - C_HALF * ( u*var_e(1,l) + v*var_e(2,l) ) );

    // Compute the normal velocity component.
    const as3double un = u * nx + v * ny;

    // Subtract the normal components to obtain the tangential ones.
    u -= un * nx;
    v -= un * ny;
   
    // Assemble the boundary state.
    var_b(0,l) = var_e(0,l);
    var_b(1,l) = var_e(0,l) * u;
    var_b(2,l) = var_e(0,l) * v;
    var_b(3,l) = p * C_OVGM1 
               + C_HALF * var_e(0,l) * ( u*u + v*v );
  }
}
