#pragma once

/*!
 * @file
 * @brief Paired-zone surface interpolation and conservative interface residuals.
 */


#include <cstddef>
#include <memory>
#include <string>

#include "option_structure.hpp"
#include "factory_structure.hpp"
#include "config_structure.hpp"
#include "geometry_structure.hpp"
#include "solver_structure.hpp"
#include "riemann_solver_structure.hpp"
#include <functional> 

// Forward declaration to avoid compiler issues.
class CMultizoneSolver;
class ISolver;


/*!
 * @brief An interface class used for the zone interface specification.
 */
class IInterface
{
	public:

    /*!
     * @brief Constructor that initializes IInterface from the specified data.
     *
     * @param[in] config_container configuration/dictionary container.
     * @param[in] multizone_solver_container solvers and residual storage for all zones.
     * @param[in] interface_batch paired-zone interface face batch.
     */
    IInterface(const CConfig               *config_container,
               const CMultizoneSolver      *multizone_solver_container,
               const CInterfaceFaceBatch   &interface_batch);

	

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~IInterface(void);

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
    virtual void ComputeInterfaceResidual(CMultizoneSolver            *multizone_solver_container,
                                          const CInterfaceFaceBatch   &face_batch,
                                          CFlattenedFaceIndex          face_info,
                                          CPoolMatrixAS3<as3double>   &workarray,
                                          as3double                    localtime) const = 0;

		/*!
		 * @brief Getter function which returns the name of the owner marker.
		 *
		 * @return mIName.
		 */
		const std::string &GetIName(void) const {return mIName;}

		/*!
		 * @brief Getter function which returns the name of the matching marker.
		 *
		 * @return mJName.
		 */
		const std::string &GetJName(void) const {return mJName;}

		/*!
		 * @brief Getter function which returns the ID of the owner zone.
		 *
		 * @return mIZone.
		 */
		std::size_t GetIZone(void) const {return mIZone;}

		/*!
		 * @brief Getter function which returns the ID of the matching zone.
		 *
		 * @return mJZone.
		 */
		std::size_t GetJZone(void) const {return mJZone;}

		/*!
		 * @brief Getter function which returns the number of faces.
		 *
		 * @return mNFace.
		 */
		std::size_t GetnFace(void) const {return mNFace;}

	protected:	
		std::string    mIName;  ///< Name of the owner interface marker.
		std::string    mJName;  ///< Name of the matching interface marker.
		std::size_t mIZone;     ///< Zone ID of the owner interface marker.
		std::size_t mJZone;     ///< Zone ID of the matching interface marker.
		EFaceLocation  mIFace;  ///< Face ID of the owner interface marker.
		EFaceLocation  mJFace;  ///< Face ID of the matching interface marker.
		
		std::size_t         mNFace;  ///< Number of faces in this collection or interface.
		std::size_t mNInt1D;         ///< Number of integration points on this marker.	

    std::size_t mNPolyI;  ///< Solution polynomial order on the owner side.
    std::size_t mNPolyJ;  ///< Solution polynomial order on the neighbor side.

		CMatrixAS3<as3double>           mWInt1D;                   ///< Integration weights on the reference element in 1D.
		std::unique_ptr<ITensorProduct> mITensorProductContainer;  ///< Tensor-product container of the iZone.
		std::unique_ptr<ITensorProduct> mJTensorProductContainer;  ///< Tensor-product container of the jZone.
		std::unique_ptr<IRiemannSolver> mRiemannSolverContainer;   ///< Riemann solver container.
	

		/*!
		 * @brief Function which returns the selected owner-side interpolation callable.
		 *
		 * @return selected interpolation or residual callable.
		 *
		 * @note Legacy declaration without a definition in the current source tree.
		 */
		inline auto GetFuncPointerInterpFaceI(void);

		/*!
		 * @brief Function which returns the selected neighbor-side interpolation callable.
		 *
		 * @return selected interpolation or residual callable.
		 *
		 * @note Legacy declaration without a definition in the current source tree.
		 */
		inline auto GetFuncPointerInterpFaceJ(void);

		/*!
		 * @brief Function which returns the selected owner-side residual callable.
		 *
		 * @return selected interpolation or residual callable.
		 *
		 * @note Legacy declaration without a definition in the current source tree.
		 */
		inline auto GetFuncPointerResidualFaceI(void);

		/*!
		 * @brief Function which returns the selected neighbor-side residual callable.
		 *
		 * @return selected interpolation or residual callable.
		 *
		 * @note Legacy declaration without a definition in the current source tree.
		 */
		inline auto GetFuncPointerResidualFaceJ(void);

	private:

		// Disable default constructor.
		IInterface(void) = delete;
		// Disable default copy constructor.
		IInterface(const IInterface&) = delete;
		// Disable default copy operator.
		IInterface& operator=(IInterface&) = delete;
};

//-----------------------------------------------------------------------------------

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



