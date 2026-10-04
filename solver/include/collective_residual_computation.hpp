#pragma once

#include <cassert>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"

// Forward declaration to avoid compiler problems.
class CGeometry;
class CMultizoneSolver;


namespace NResidualComputation
{
  // Collective contract for every stage (when OpenMP is enabled):
  // All threads in the enclosing parallel team call in the same order,
  // outside single/master, explicit tasks, and other worksharing loops.
  // Each stage includes an implicit end-of-loop barrier.
  // Preserve the original per-thread workspace ownership arrangement.
  
  void ComputeVolumeResidualsCollective(const CGeometry            *geometry_container,
                                        CMultizoneSolver           *multizone_solver_container,
                                        CPoolMatrixAS3<as3double>  &workarray,
                                        as3double                   localtime);
  
  void ComputeIFaceResidualsCollective(const CGeometry           *geometry_container,
                                       CMultizoneSolver          *multizone_solver_container,
                                       CPoolMatrixAS3<as3double> &workarray,
                                       as3double                  localtime);
   
  void ComputeJFaceResidualsCollective(const CGeometry           *geometry_container,
                                       CMultizoneSolver          *multizone_solver_container,
                                       CPoolMatrixAS3<as3double> &workarray,
                                       as3double                  localtime);
  
  void AccumulateIFaceResidualsCollective(const CGeometry  *geometry_container,
                                          CMultizoneSolver *multizone_solver_container);
  
  void AccumulateJFaceResidualsCollective(const CGeometry  *geometry_container,
                                          CMultizoneSolver *multizone_solver_container);
  
  void ApplyInverseMassMatricesCollective(const CGeometry           *geometry_container,
                                          CMultizoneSolver          *multizone_solver_container,
                                          CPoolMatrixAS3<as3double> &workarray);

  inline void AssertParallelRegion(void)
  {
#if defined(HAVE_OPENMP) && !defined(NDEBUG)
    assert( omp_get_level() > 0 );
#endif
  }
}
