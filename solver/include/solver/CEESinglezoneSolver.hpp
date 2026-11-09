#pragma once

/*!
 * @file
 * @brief Euler equations solver for one grid zone.
 */


#include "solver/ISinglezoneSolver.hpp"
#include "data/CPoolMatrixAS3.hpp"


/*!
 * @brief A class for a solver specification based on the (non-linear) Euler equations.
 */
class CEESinglezoneSolver final: public ISinglezoneSolver
{
	public:

		/*!
		 * @brief Constructor that initializes CEESinglezoneSolver from the specified data.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 */
		CEESinglezoneSolver(const CConfig             *config_container,
				                const CSinglezoneGeometry *singlezone_geometry_container);
		

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CEESinglezoneSolver(void) override;

		/*!
		 * @brief Getter function which returns the governing-equation type.
		 *
		 * @return implemented governing-equation type.
		 */
		ETypeSolver GetTypeSolver(void) const override {return ETypeSolver::EE;}

		/*!
		 * @brief Function that initializes the physical elements.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] multizone_geometry_container input geometry container.
		 */
		void InitPhysicalElements(const CConfig            *config_container,
				                      const CMultizoneGeometry *multizone_geometry_container) override;

		/*!
		 * @brief Function that computes the element-local volume and physical face terms for the EE.
		 *
		 * @param[in] singlezone_geometry_container geometry of the current grid zone.
		 * @param[in,out] workarray memory for the working array.
		 * @param[in] localtime current physical time.
		 * @param[in] iElem current element index.
		 *
		 * @note Overwrites the local residual. Weak form includes only volume terms; strong
		 *   form also includes the four physical face fluxes. Coupled gathering and inverse mass
		 *   application are performed by the element residual operation.
		 */
		void ComputeElementLocalResidual(CSinglezoneGeometry       *singlezone_geometry_container,
				                       CPoolMatrixAS3<as3double> &workarray,
															 as3double                  localtime,
															 std::size_t                     iElem) override;


		/*!
		 * @brief Function that computes one same-zone I-face flux and writes the two adjacent face
		 *   residual buffers.
		 *
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in,out] workarray workspace pool owned by the calling thread.
		 * @param[in] localtime physical time at which the residual is evaluated.
		 * @param[in] iElemL zone-local minus-side element adjacent to the I-face.
		 */
		void ComputeElementCoupledResidualIDir(const CSinglezoneGeometry *singlezone_geometry_container,
		                                CPoolMatrixAS3<as3double> &workarray,
																	  as3double                  localtime,
																	  std::size_t                     iElemL) override final;

		/*!
		 * @brief Function that computes one same-zone J-face flux and writes the two adjacent face
		 *   residual buffers.
		 *
		 * @param[in] singlezone_geometry_container geometry of the current zone.
		 * @param[in,out] workarray workspace pool owned by the calling thread.
		 * @param[in] localtime physical time at which the residual is evaluated.
		 * @param[in] iElemB zone-local minus-side element adjacent to the J-face.
		 */
		void ComputeElementCoupledResidualJDir(const CSinglezoneGeometry *singlezone_geometry_container,
		                                CPoolMatrixAS3<as3double> &workarray,
																	  as3double                  localtime,
																	  std::size_t                     iElemB) override final;

		/*!
		 * @brief Getter function which returns the number of solution variables.
		 *
		 * @return mNVar.
		 */
		std::size_t GetnVar(void) const override {return mNVar;}

	protected:

	private:
		inline static constexpr std::size_t mNVar = 4; ///< Number of working variables

};
