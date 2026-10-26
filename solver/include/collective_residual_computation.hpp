#pragma once

/*!
 * @file
 * @brief Collective OpenMP surface sweeps and element-local residual assembly.
 */


#include <cstddef>

#include <cassert>

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"

// Forward declaration to avoid compiler problems.
class CMultizoneGeometry;
class CMultizoneSolver;
class ISolver;


/*!
 * @brief Surface and element residual assembly with explicit OpenMP synchronization contracts.
 *
 * @note With OpenMP enabled, collective calls are encountered by every thread in the same enclosing team, in
 *   the same order.
 * Call outside single/master regions, explicit tasks and other worksharing loops. Workspace pools are private
 *   to each thread.
 */
namespace NResidualComputation
{

  /*!
   * @brief Function that computes both directional surface residual sweeps without an end barrier.
   *
   * @param[in] multizone_geometry_container geometry of all grid zones.
   * @param[in,out] multizone_solver_container solvers and residual storage for all zones.
   * @param[in,out] workarray workspace pool owned by the calling thread.
   * @param[in] localtime physical time at which the residual is evaluated.
   *
   * @note All threads in the enclosing OpenMP team must call collectively in the same order.
   *
   * @note Both directional worksharing loops use nowait. Synchronize before consuming the face buffers or
   *   updating the solution.
   */
  
  void ComputeSurfaceResidualsCollectiveNowait(const CMultizoneGeometry  *multizone_geometry_container,
                                               CMultizoneSolver          *multizone_solver_container,
                                               CPoolMatrixAS3<as3double> &workarray,
                                               as3double                  localtime);
  

  /*!
   * @brief Function that collectively computes volume residuals, gathers surface contributions and applies
   *   inverse mass matrices.
   *
   * @param[in] multizone_geometry_container geometry of all grid zones.
   * @param[in,out] multizone_solver_container solvers and residual storage for all zones.
   * @param[in,out] workarray workspace pool owned by the calling thread.
   * @param[in] localtime physical time at which the residual is evaluated.
   *
   * @note All threads in the enclosing OpenMP team must call collectively in the same order.
   *
   * @note An entry barrier waits for surface writes. The element worksharing loop has an implicit exit barrier.
   */
  void ComputeElementResidualsCollective(const CMultizoneGeometry  *multizone_geometry_container,
                                         CMultizoneSolver          *multizone_solver_container,
                                         CPoolMatrixAS3<as3double> &workarray,
                                         as3double                  localtime);

  /*!
   * @brief Function that computes one element volume residual, gathers its four surface residuals and applies
   *   its inverse mass matrix.
   *
   * @param[in] singlezone_geometry_container geometry of the current zone.
   * @param[in,out] singlezone_solver_container solver and element storage for the current zone.
   * @param[in,out] workarray workspace pool owned by the calling thread.
   * @param[in] localtime physical time at which the residual is evaluated.
   * @param[in] iElem zone-local element index.
   *
   * @note No worksharing or barrier is performed. Complete surface writes before calling; assign each element
   *   to one thread.
   */
  void ComputeElementResidual(CSinglezoneGeometry       *singlezone_geometry_container,
                               ISolver                   *singlezone_solver_container,
                               CPoolMatrixAS3<as3double> &workarray,
                               as3double                  localtime,
                               std::size_t                     iElem);

  /*!
   * @brief Function that checks collective entry from an OpenMP parallel region in assertion-enabled builds.
   */
  inline void AssertParallelRegion(void)
  {
#if defined(HAVE_OPENMP) && !defined(NDEBUG)
    assert( omp_get_level() > 0 );
#endif
  }
}
