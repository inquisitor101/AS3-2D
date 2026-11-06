#pragma once

/*!
 * @file
 * @brief Euler fluxes on paired zone edges.
 */


#include <functional>

#include "interface/IInterface.hpp"
#include "data/CPoolMatrixAS3.hpp"


/*!
 * @brief A class for an interface specification based on the (non-linear) Euler equations.
 */
class CEEInterface final: public IInterface
{
  private:

    /*!
     * @brief Callable accepting a solution and optional value, r-derivative and s-derivative outputs.
     */
    using AInterpolateSurface = std::function<void(const as3double*,
                                                         as3double*,
                                                         as3double*,
                                                         as3double*)>;
	

    /*!
     * @brief Callable that lifts weighted value and derivative inputs into a face residual.
     */
    using AComputeResidualFace = std::function<void(const as3double*,
                                                          as3double*,
                                                          as3double*,
                                                          as3double*)>;
  public:

    /*!
     * @brief Constructor that initializes CEEInterface from the specified data.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] multizone_solver_container solvers and residual storage for all zones.
     * @param[in] interface_batch paired-zone interface face batch.
     */
    CEEInterface(const CConfig               *config_container,
                 const CMultizoneSolver      *multizone_solver_container,
                 const CInterfaceFaceBatch   &interface_batch);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CEEInterface(void) override;

    /*!
     * @brief Function that computes one paired-face flux and writes each side to its own local surface residual
     *   buffer.
     *
     * @param[in,out] multizone_solver_container solvers and residual storage for all zones.
     * @param[in] face_batch connectivity and orientation of the current face batch.
     * @param[in] face_info face category and category-local batch and face indices.
     * @param[in,out] workarray workspace pool owned by the calling thread.
     * @param[in] localtime physical time at which the residual is evaluated.
     */
    void ComputeInterfaceResidual(CMultizoneSolver            *multizone_solver_container,
                                  const CInterfaceFaceBatch   &face_batch,
                                  CFlattenedFaceIndex          face_info,
                                  CPoolMatrixAS3<as3double>   &workarray,
                                  as3double                    localtime) const override final;
	protected:

	private:
    std::size_t mNVar; ///< Number of solution variables.

    AComputeResidualFace mFComputeResidualFaceI;  ///< Selected owner-side residual lifting kernel.
    AComputeResidualFace mFComputeResidualFaceJ;  ///< Selected neighbor-side residual lifting kernel.

    AInterpolateSurface  mFInterpolateSurfaceI;  ///< Selected owner-side trace interpolation kernel.
    AInterpolateSurface  mFInterpolateSurfaceJ;  ///< Selected neighbor-side trace interpolation kernel.

    /*!
     * @brief Function that selects face interpolation and residual lifting kernels for the local edges.
     */
    void InitializeComputeKernels(void);
};
