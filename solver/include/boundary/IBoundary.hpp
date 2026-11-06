#pragma once

/*!
 * @file
 * @brief Physical boundary interface and common residual evaluation.
 */


#include <functional>
#include <memory>

#include "common/AS3Options.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CPoolMatrixAS3.hpp"

class CConfig;
class CSinglezoneGeometry;
class ISinglezoneSolver;
class CBoundaryFaceBatch;
class IRiemannSolver;
struct CFlattenedFaceIndex;
struct IBoundaryParamMarker;


/*!
 * @brief An interface class used for the boundary condition.
 */
class IBoundary
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
     * @brief Constructor that initializes IBoundary from the specified data.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] singlezone_geometry_container geometry of the current zone.
     * @param[in] singlezone_solver_container solver and element storage for the current zone.
     * @param[in] boundary_batch physical boundary face batch.
     * @param[in] param_bc parameters associated with this physical boundary marker.
     */
    IBoundary(const CConfig              *config_container,
              const CSinglezoneGeometry  *singlezone_geometry_container,
              const ISinglezoneSolver    *singlezone_solver_container,
              const CBoundaryFaceBatch   &boundary_batch,
              const IBoundaryParamMarker *param_bc);

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~IBoundary(void);


		/*!
		 * @brief Getter function which returns the boundary-condition type implemented by this object.
		 *
		 * @return boundary-condition type implemented by this object.
		 */
		virtual ETypeBC GetTypeBC(void) const = 0;

    /*!
     * @brief Function that computes one boundary flux and writes the adjacent element face residual buffer.
     *
     * @param[in] singlezone_geometry_container geometry of the current zone.
     * @param[in,out] singlezone_solver_container solver and element storage for the current zone.
     * @param[in] face_batch connectivity and orientation of the current face batch.
     * @param[in] face_info face category and category-local batch and face indices.
     * @param[in,out] workarray workspace pool owned by the calling thread.
     * @param[in] localtime physical time at which the residual is evaluated.
     */
    void ComputeBoundaryResidual(const CSinglezoneGeometry *singlezone_geometry_container,
                                 ISinglezoneSolver         *singlezone_solver_container,
                                 const CBoundaryFaceBatch  &face_batch,
                                 CFlattenedFaceIndex        face_info,
                                 CPoolMatrixAS3<as3double> &workarray,
                                 as3double                  localtime) const;

	protected:
		
    std::unique_ptr<IRiemannSolver> mRiemannSolverContainer;  ///< Riemann solver container.

    AInterpolateSurface  mFInterpolateSurface;   ///< Selected interpolation kernel for the local boundary edge.
    AComputeResidualFace mFComputeResidualFace;  ///< Selected kernel for lifting this face flux into an element residual.

    /*!
     * @brief Function that constructs the exterior conservative state used by the boundary Riemann solver.
     *
     * @param[in] met surface metric rows: length factor, outward normal x, outward normal y.
     * @param[in] var_e interior conservative state at the boundary integration points.
     * @param[out] var_b exterior conservative state passed to the Riemann solver.
     */
    virtual void ComputeBoundaryState(const CMatrixAS3<as3double>     &met,
                                      const CWorkMatrixAS3<as3double> &var_e,
                                      CWorkMatrixAS3<as3double>       &var_b) const = 0;



	private:
		// Disable default constructor.
		IBoundary(void) = delete;
		// Disable default copy constructor.
		IBoundary(const IBoundary&) = delete;
		// Disable default copy operator.
		IBoundary& operator=(IBoundary&) = delete;


    /*!
     * @brief Function that selects face interpolation and residual lifting kernels for the local edges.
     *
     * @param[in,out] singlezone_solver_container solver and element storage for the current zone.
     * @param[in] face_location local edge: IMIN, IMAX, JMIN or JMAX.
     */
    void InitializeComputeKernels(const ISinglezoneSolver *singlezone_solver_container,
                                  EFaceLocation            face_location);
};
