#pragma once

/*!
 * @file
 * @brief Riemann inlet with prescribed static mass flow.
 */


#include <array>

#include "boundary/IBoundary.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CWorkMatrixAS3.hpp"


/*!
 * @brief An inlet that supplies prescribed density and velocity to the Riemann solver.
 */
class CRiemannStaticMassInletBoundary final : public IBoundary
{
	public:
	

		/*!
		 * @brief Constructor that initializes CRiemannStaticMassInletBoundary from the specified data.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in] singlezone_solver_container solver and element storage for the current zone.
		 * @param[in] boundary_batch physical boundary face batch.
		 * @param[in] param_bc parameters associated with this physical boundary marker.
		 */
		CRiemannStaticMassInletBoundary(const CConfig              *config_container,
                                    const CSinglezoneGeometry  *singlezone_geometry_container,
                                    const ISinglezoneSolver    *singlezone_solver_container,
                                    const CBoundaryFaceBatch   &boundary_batch,
                                    const IBoundaryParamMarker *param_bc);

    /*!
     * @brief Destructor that releases the resources owned by this object.
     */
    ~CRiemannStaticMassInletBoundary(void) override;

    /*!
     * @brief Getter function which returns the boundary-condition type implemented by this marker or boundary.
     *
     * @return boundary-condition type implemented by this object.
     */
    ETypeBC GetTypeBC(void) const override { return ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN; }

  protected:

  private:
    as3double                mPrescribedDensity;   ///< Prescribed inlet static density.
    std::array<as3double, 2> mPrescribedVelocity;  ///< Prescribed inlet Cartesian velocity components.

    /*!
     * @brief Function that constructs the exterior conservative state used by the boundary Riemann solver.
     *
     * @param[in] met surface metric rows: length factor, outward normal x, outward normal y.
     * @param[in] var_e interior conservative state at the boundary integration points.
     * @param[out] var_b exterior conservative state passed to the Riemann solver.
     */
    void ComputeBoundaryState(const CMatrixAS3<as3double>     &met,
                              const CWorkMatrixAS3<as3double> &var_e,
                              CWorkMatrixAS3<as3double>       &var_b) const override;
};
