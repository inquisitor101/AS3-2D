#pragma once

/*!
 * @file
 * @brief Interface for residual evaluation between two grid zones.
 */


#include <cstddef>
#include <memory>
#include <string>

#include "common/AS3Options.hpp"
#include "data/CMatrixAS3.hpp"
#include "data/CPoolMatrixAS3.hpp"

class CConfig;
class CInterfaceFaceBatch;
class CMultizoneGeometry;
class CMultizoneSolver;
class ITensorProduct;
class IRiemannSolver;
struct CFlattenedFaceIndex;


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
