#pragma once

/*!
 * @file
 * @brief Singlezone solver interface and element storage.
 */


#include <cstddef>
#include <memory>

#include "common/AS3Options.hpp"
#include "data/CPoolMatrixAS3.hpp"

class CConfig;
class CMultizoneGeometry;
class CSinglezoneGeometry;
class ITensorProduct;
class IRiemannSolver;
class CStandardElement;
class CPhysicalElement;


/*!
 * @brief Interface and owned solution data for the solver in one grid zone.
 */
class ISinglezoneSolver
{
	public:
		

		/*!
		 * @brief Constructor that initializes ISinglezoneSolver from the specified data.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in] nVar number of solution variables.
		 */
		ISinglezoneSolver(const CConfig             *config_container,
				              const CSinglezoneGeometry *singlezone_geometry_container,
                      std::size_t                nVar);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		virtual ~ISinglezoneSolver(void);

		/*!
		 * @brief Getter function which returns the governing-equation type.
		 *
		 * @return implemented governing-equation type.
		 */
		virtual ETypeSolver GetTypeSolver(void) const = 0;

		/*!
		 * @brief Pure virtual function that initializes the physical elements. Must be implemented in a
		 *   derived class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input multizone geometry container.
		 */
		virtual void InitPhysicalElements(const CConfig            *config_container,
				                              const CMultizoneGeometry *multizone_geometry_container) = 0;

		/*!
		 * @brief Pure virtual function that computes the element-local volume and physical face terms.
		 *
		 * @param[in] singlezone_container geometry of the current grid zone.
		 * @param[in,out] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 * @param[in] iElem current element index.
		 *
		 * @note Overwrites the local residual. Weak form includes only volume terms; strong form also
		 *   includes the four physical face fluxes. Coupled gathering and inverse mass application are
		 *   performed by the element residual operation.
		 */
		virtual void ComputeElementLocalResidual(CSinglezoneGeometry       *singlezone_container,
																			 CPoolMatrixAS3<as3double> &workarray,
																			 as3double                  localtime,
																			 std::size_t                     iElem) = 0;

		/*!
		 * @brief Function that computes one same-zone I-face flux and writes the two adjacent face
		 *   residual buffers.
		 *
		 * @param[in] singlezone_container geometry of the current zone.
		 * @param[in,out] workarray workspace pool owned by the calling thread.
		 * @param[in] localtime physical time at which the residual is evaluated.
		 * @param[in] iElemL zone-local minus-side element adjacent to the I-face.
		 */
		virtual void ComputeElementCoupledResidualIDir(const CSinglezoneGeometry *singlezone_container,
				                                    CPoolMatrixAS3<as3double> &workarray,
																					  as3double                  localtime,
																					  std::size_t                     iElemL) = 0;
		

    /*!
     * @brief Function that computes one same-zone J-face flux and writes the two adjacent face residual
     *   buffers.
     *
     * @param[in] singlezone_container geometry of the current zone.
     * @param[in,out] workarray workspace pool owned by the calling thread.
     * @param[in] localtime physical time at which the residual is evaluated.
     * @param[in] iElemB zone-local minus-side element adjacent to the J-face.
     */
    virtual void ComputeElementCoupledResidualJDir(const CSinglezoneGeometry *singlezone_container,
				                                    CPoolMatrixAS3<as3double> &workarray,
																					  as3double                  localtime,
																					  std::size_t                     iElemB) = 0;


		/*!
		 * @brief Getter function which returns the number of solution variables.
		 *
		 * @return number of solution variables.
		 */
		virtual std::size_t GetnVar(void) const = 0;

		/*!
		 * @brief Getter function which returns the value of mZoneID.
		 *
		 * @return mZoneID.
		 */
		std::size_t GetZoneID(void) const {return mZoneID;}

		/*!
		 * @brief Getter function which returns the tensor-product container of this zone.
		 *
		 * @return mTensorProductContainer.
		 */
		ITensorProduct *GetTensorProduct(void) const {return mTensorProductContainer.get();}

		/*!
		 * @brief Getter function which returns the Riemann solver container of this zone.
		 *
		 * @return mRiemannSolverContainer.
		 */
		IRiemannSolver *GetRiemannSolver(void) const {return mRiemannSolverContainer.get();}

		/*!
		 * @brief Getter function which returns the standard element container of this zone.
		 *
		 * @return mStandardElementContainer.
		 */
		const CStandardElement *GetStandardElement(void) const {return mStandardElementContainer.get();}

		/*!
		 * @brief Getter function which returns the entire physical element container.
		 *
		 * @return mPhysicalElementContainer.
		 */
		as3vector1d<std::unique_ptr<CPhysicalElement>> &GetPhysicalElement(void) {return mPhysicalElementContainer;}

		/*!
		 * @brief Getter function which returns a specific physical element container.
		 *
		 * @param[in] index index of the physical element.
		 *
		 * @return mPhysicalElementContainer[index].
		 */
		CPhysicalElement *GetPhysicalElement(std::size_t index) const {return mPhysicalElementContainer[index].get();}


	protected:
		const std::size_t                              mZoneID;                    ///< Zone ID of this container.
		const ETypeDGFormulation                       mTypeDGFormulation;         ///< Weak or strong local residual operator.
		std::unique_ptr<ITensorProduct>                mTensorProductContainer;    ///< Tensor product container.
		std::unique_ptr<IRiemannSolver>                mRiemannSolverContainer;    ///< Riemann solver container.
		std::unique_ptr<CStandardElement>              mStandardElementContainer;  ///< Standard element container.	
		as3vector1d<std::unique_ptr<CPhysicalElement>> mPhysicalElementContainer;  ///< Physical element container.
	
	private:
		// Disable default constructor.
		ISinglezoneSolver(void) = delete;
		// Disable default copy constructor.
		ISinglezoneSolver(const ISinglezoneSolver&) = delete;
		// Disable default copy operator.
		ISinglezoneSolver& operator=(ISinglezoneSolver&) = delete;
};
