#pragma once 

#include "option_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "marker_structure.hpp"
#include "factory_structure.hpp"
#include "riemann_solver_structure.hpp"
#include "solver_structure.hpp"

// Forward declaration to avoid compiler issues.
class ISolver;



/*!
 * @brief An interface class used for the boundary condition.
 */
class IBoundary
{
  private:

    using AInterpolateSurface = std::function<void(const as3double*,
                                                         as3double*,
                                                         as3double*,
                                                         as3double*)>;

    using AComputeResidualFace = std::function<void(const as3double*,
                                                          as3double*,
                                                          as3double*,
                                                          as3double*)>;
	public:
	
    IBoundary(const CConfig              *config_container,
				      const CSinglezoneGeometry  *singlezone_geometry_container,
              const ISolver              *singlezone_solver_container,
              const CBoundaryFacesFamily &boundary_family,
              const IBoundaryParamMarker *param_bc);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~IBoundary(void);


		/*!
		 * @brief Pure virtual function that returns the type of BC. Must be overridden.
		 */
		virtual ETypeBC GetTypeBC(void) const = 0;

    void ComputeBoundaryResidual(const CSinglezoneGeometry  *singlezone_geometry_container,
                                 ISolver                    *singlezone_solver_container,
                                 const CBoundaryFacesFamily &family_face,
                                 CFlattenedFaceIndex         face_info,
                                 CPoolMatrixAS3<as3double>  &workarray,
                                 as3double                   localtime) const;

	protected:
		
    std::unique_ptr<IRiemannSolver> mRiemannSolverContainer;  ///< Riemann solver container.

    AInterpolateSurface  mFInterpolateSurface;
    AComputeResidualFace mFComputeResidualFace;

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


    void InitializeComputeKernels(const ISolver *singlezone_solver_container,
                                  EFaceLocation  face_location);
};

//-----------------------------------------------------------------------------------

class CSymmetryBoundary : public IBoundary
{
	public:
	
		CSymmetryBoundary(const CConfig              *config_container,
				              const CSinglezoneGeometry  *singlezone_geometry_container,
                      const ISolver              *singlezone_solver_container,
                      const CBoundaryFacesFamily &boundary_family,
                      const IBoundaryParamMarker *param_bc);

    ~CSymmetryBoundary(void) final;

    ETypeBC GetTypeBC(void) const final { return ETypeBC::SYMMETRY; }

  protected:

  private:

    void ComputeBoundaryState(const CMatrixAS3<as3double>     &met,
                              const CWorkMatrixAS3<as3double> &var_e,
                              CWorkMatrixAS3<as3double>       &var_b) const;

};


//-----------------------------------------------------------------------------------

class CRiemannStaticMassInletBoundary : public IBoundary
{
	public:
	
		CRiemannStaticMassInletBoundary(const CConfig              *config_container,
				                            const CSinglezoneGeometry  *singlezone_geometry_container,
                                    const ISolver              *singlezone_solver_container,
                                    const CBoundaryFacesFamily &boundary_family,
                                    const IBoundaryParamMarker *param_bc);

    ~CRiemannStaticMassInletBoundary(void) final;

    ETypeBC GetTypeBC(void) const final { return ETypeBC::STATIC_MASS_FLOW_INLET_RIEMANN; }

  protected:

  private:
    as3double                mPrescribedDensity;
    std::array<as3double, 2> mPrescribedVelocity;

    void ComputeBoundaryState(const CMatrixAS3<as3double>     &met,
                              const CWorkMatrixAS3<as3double> &var_e,
                              CWorkMatrixAS3<as3double>       &var_b) const;
};


//-----------------------------------------------------------------------------------

class CRiemannPressureOutletBoundary : public IBoundary
{
	public:
	
		CRiemannPressureOutletBoundary(const CConfig              *config_container,
				                           const CSinglezoneGeometry  *singlezone_geometry_container,
                                   const ISolver              *singlezone_solver_container,
                                   const CBoundaryFacesFamily &boundary_family,
                                   const IBoundaryParamMarker *param_bc);

    ~CRiemannPressureOutletBoundary(void) final;

    ETypeBC GetTypeBC(void) const final { return ETypeBC::PRESSURE_OUTLET_RIEMANN; }

  protected:

  private:
    as3double mPrescribedPressure;

    void ComputeBoundaryState(const CMatrixAS3<as3double>     &met,
                              const CWorkMatrixAS3<as3double> &var_e,
                              CWorkMatrixAS3<as3double>       &var_b) const;

};  

