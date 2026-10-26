#pragma once

/*!
 * @file
 * @brief Runtime and compile-time tensor interpolation, residual lifting and mass application.
 */


#include <cstddef>

#include "option_structure.hpp"
#include "standard_element_structure.hpp"
#include "compile_time_dimensions.hpp"



/*!
 * @brief An interface for runtime tensor products and compile-time tensor specializations.
 */
class ITensorProduct
{
	public:

		/*!
		 * @brief Constructor of ITensorProduct, which serves as an interface for the tensor product.
		 */
		ITensorProduct(void) = default;
		

		/*!
		 * @brief Virtual destructor for destruction through the component interface.
		 */
		virtual ~ITensorProduct(void) = default;

    /*!
     * @brief Getter function which returns the number of solution points in one reference direction.
     *
     * @return number of solution points in one reference direction.
     */
    virtual std::size_t GetnSol1D(void) const = 0;

    /*!
     * @brief Getter function which returns the number of integration points in one reference direction.
     *
     * @return number of integration points in one reference direction.
     */
    virtual std::size_t GetnInt1D(void) const = 0;

    /*!
     * @brief Getter function which returns the number of solution variables.
     *
     * @return number of solution variables.
     */
    virtual std::size_t GetnVar(void)   const = 0;


		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   volume integration points.
		 *
		 * @param[in] K number of solution points in 1D.
		 * @param[in] N number of variables.
		 * @param[in] M number of integration points in 1D.
		 * @param[in] A transposed interpolation operator [K][M].
		 * @param[in] ADer transposed parametric differentiation operator [K][M].
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M*M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M*M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M*M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void RuntimeVolume(std::size_t      K,
				               std::size_t      N,
 				               std::size_t      M,
 				               const as3double *A,
 				               const as3double *ADer,
 				               const as3double *B,
 				               as3double       *C,
 				               as3double       *CDerR,
 				               as3double       *CDerS);

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMIN surface integration points.
		 *
		 * @param[in] K number of solution points in 1D.
		 * @param[in] N number of variables.
		 * @param[in] M number of integration points in 1D.
		 * @param[in] A transposed interpolation operator [K][M].
		 * @param[in] ADer transposed parametric differentiation operator [K][M].
		 * @param[in] ADerFace reference-coordinate derivative at the selected endpoint, containing K
		 *   entries.
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void RuntimeSurfaceIMIN(std::size_t      K,
				                    std::size_t      N,
 				                    std::size_t      M,
 				                    const as3double *A,
 				                    const as3double *ADer,
													  const as3double *ADerFace,
 				                    const as3double *B,
 				                    as3double       *C,
 				                    as3double       *CDerR,
 				                    as3double       *CDerS);

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMAX surface integration points.
		 *
		 * @param[in] K number of solution points in 1D.
		 * @param[in] N number of variables.
		 * @param[in] M number of integration points in 1D.
		 * @param[in] A transposed interpolation operator [K][M].
		 * @param[in] ADer transposed parametric differentiation operator [K][M].
		 * @param[in] ADerFace reference-coordinate derivative at the selected endpoint, containing K
		 *   entries.
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void RuntimeSurfaceIMAX(std::size_t      K,
				                    std::size_t      N,
 				                    std::size_t      M,
 				                    const as3double *A,
 				                    const as3double *ADer,
													  const as3double *ADerFace,
 				                    const as3double *B,
 				                    as3double       *C,
 				                    as3double       *CDerR,
 				                    as3double       *CDerS);

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMIN surface integration points.
		 *
		 * @param[in] K number of solution points in 1D.
		 * @param[in] N number of variables.
		 * @param[in] M number of integration points in 1D.
		 * @param[in] A transposed interpolation operator [K][M].
		 * @param[in] ADer transposed parametric differentiation operator [K][M].
		 * @param[in] ADerFace reference-coordinate derivative at the selected endpoint, containing K
		 *   entries.
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void RuntimeSurfaceJMIN(std::size_t      K,
				                    std::size_t      N,
 				                    std::size_t      M,
 				                    const as3double *A,
 				                    const as3double *ADer,
													  const as3double *ADerFace,
 				                    const as3double *B,
 				                    as3double       *C,
 				                    as3double       *CDerR,
 				                    as3double       *CDerS);

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMAX surface integration points.
		 *
		 * @param[in] K number of solution points in 1D.
		 * @param[in] N number of variables.
		 * @param[in] M number of integration points in 1D.
		 * @param[in] A transposed interpolation operator [K][M].
		 * @param[in] ADer transposed parametric differentiation operator [K][M].
		 * @param[in] ADerFace reference-coordinate derivative at the selected endpoint, containing K
		 *   entries.
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void RuntimeSurfaceJMAX(std::size_t      K,
				                    std::size_t      N,
 				                    std::size_t      M,
 				                    const as3double *A,
 				                    const as3double *ADer,
													  const as3double *ADerFace,
 				                    const as3double *B,
 				                    as3double       *C,
 				                    as3double       *CDerR,
 				                    as3double       *CDerS);




		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   volume integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M*M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M*M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M*M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void PartialCompileTimeVolume(std::size_t      n,
                                          const as3double *B,
								                          as3double       *C,
								                          as3double       *CDerR,
								                          as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMIN surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void PartialCompileTimeSurfaceIMIN(std::size_t      n,
                                               const as3double *B,
								                               as3double       *C,
								                               as3double       *CDerR,
								                               as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMAX surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void PartialCompileTimeSurfaceIMAX(std::size_t      n,
                                               const as3double *B,
								                               as3double       *C,
								                               as3double       *CDerR,
								                               as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMIN surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void PartialCompileTimeSurfaceJMIN(std::size_t      n,
                                               const as3double *B,
								                               as3double       *C,
								                               as3double       *CDerR,
								                               as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMAX surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void PartialCompileTimeSurfaceJMAX(std::size_t      n,
                                               const as3double *B,
								                               as3double       *C,
								                               as3double       *CDerR,
								                               as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   volume integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M*M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M*M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M*M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeVolume(const as3double *B,
								                   as3double       *C,
								                   as3double       *CDerR,
								                   as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMIN surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeSurfaceIMIN(const as3double *B,
								                        as3double       *C,
								                        as3double       *CDerR,
								                        as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMAX surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeSurfaceIMAX(const as3double *B,
								                        as3double       *C,
								                        as3double       *CDerR,
								                        as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMIN surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeSurfaceJMIN(const as3double *B,
								                        as3double       *C,
								                        as3double       *CDerR,
								                        as3double       *CDerS) = 0;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMAX surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeSurfaceJMAX(const as3double *B,
								                        as3double       *C,
								                        as3double       *CDerR,
								                        as3double       *CDerS) = 0;

		/*!
		 * @brief Function that lifts weighted inputs from the volume integration points to the solution
		 *   DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M*M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M*M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M*M]; nullptr omits this
		 *   contribution.
		 * @param[out] C solution residual [N][K*K], overwritten by the volume contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeResidualVolume(const as3double *B,
																           const as3double *BDerR,
																           const as3double *BDerS,
																           as3double       *C) = 0;

		/*!
		 * @brief Function that lifts weighted inputs from the IMIN surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeResidualSurfaceIMIN(const as3double *B,
																                const as3double *BDerR,
																                const as3double *BDerS,
																                as3double       *C) = 0;

		/*!
		 * @brief Function that lifts weighted inputs from the IMAX surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeResidualSurfaceIMAX(const as3double *B,
																                const as3double *BDerR,
																                const as3double *BDerS,
																                as3double       *C) = 0;

		/*!
		 * @brief Function that lifts weighted inputs from the JMIN surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeResidualSurfaceJMIN(const as3double *B,
																                const as3double *BDerR,
																                const as3double *BDerS,
																                as3double       *C) = 0;

		/*!
		 * @brief Function that lifts weighted inputs from the JMAX surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		virtual void CompileTimeResidualSurfaceJMAX(const as3double *B,
                                                const as3double *BDerR,
                                                const as3double *BDerS,
                                                as3double       *C) = 0;


    /*!
     * @brief Function that replaces the element residual by its inverse-mass product.
     *
     * @param[in] minv row-major inverse mass matrix of size [K*K][K*K].
     * @param[in,out] res element residual [N][K*K], replaced by its inverse-mass product.
     * @param[out] tmp scratch storage for at least K*K entries, separate from the residual.
     *
     * @note The scratch buffer must not overlap the residual or inverse mass matrix.
     */
    virtual void CompileTimeApplyInverseMassMatrix(const as3double *minv,
                                                   as3double       *res,
                                                   as3double       *tmp) = 0;
	protected:

	private:
		inline static constexpr std::size_t mRuntimeMaxK = 10;  ///< Maximum solution-point count supported by the runtime tensor workspace.
    inline static constexpr std::size_t mRuntimeMaxM = 10;  ///< Maximum integration-point count supported by the runtime tensor workspace.
    
    // Disable default copy constructor.
		ITensorProduct(const ITensorProduct&) = delete;
		// Disable default copy operator.
		ITensorProduct& operator=(ITensorProduct&) = delete;	
};

//-----------------------------------------------------------------------------------

/*!
 * @brief Tensor-product operations specialized for K solution points, M integration points and N variables.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 */
template<std::size_t K, std::size_t M, std::size_t N>
class CTensorProduct final: public ITensorProduct
{
  private:

    /*!
     * @brief Compile-time dimensions associated with this tensor specialization.
     */
    using CDimensions = CCompileTimeDimensions<K,M,N>;
	
  public:
		

		/*!
		 * @brief Constructor which initializes this specialized tensor-product class.
		 *
		 * @param[in] standard_element pointer to the standard element of this specialization.
		 */
		explicit CTensorProduct(CStandardElement *standard_element) :
				mLagrangeInt1D        ( standard_element->GetLagrangeInt1D()         ),
				mLagrangeInt1DTrans   ( standard_element->GetLagrangeInt1DTrans()    ),
				mDerLagrangeInt1D     ( standard_element->GetDerLagrangeInt1D()      ),
				mDerLagrangeInt1DTrans( standard_element->GetDerLagrangeInt1DTrans() ),
				mDerLagrangeMinFace1D ( standard_element->GetDerLagrangeMinFace1D()  ),
				mDerLagrangeMaxFace1D ( standard_element->GetDerLagrangeMaxFace1D()  ) 
		{
			if( K != standard_element->GetnSol1D() ) ERROR("Mismatch in tensor K dimension.");
			if( M != standard_element->GetnInt1D() ) ERROR("Mismatch in tensor M dimension.");
		} 

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CTensorProduct(void) override = default;

    /*!
     * @brief Getter function which returns the number of solution points in one reference direction.
     *
     * @return number of solution points in one reference direction.
     */
    std::size_t GetnSol1D(void) const override { return CDimensions::nSol1D; }

    /*!
     * @brief Getter function which returns the number of integration points in one reference direction.
     *
     * @return number of integration points in one reference direction.
     */
    std::size_t GetnInt1D(void) const override { return CDimensions::nInt1D; }

    /*!
     * @brief Getter function which returns the number of solution variables.
     *
     * @return number of solution variables.
     */
    std::size_t GetnVar(void)   const override { return CDimensions::nVar;   }


		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   volume integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M*M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M*M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M*M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void PartialCompileTimeVolume(std::size_t      n,
                                  const as3double *B,
				                          as3double       *C,
				                          as3double       *CDerR,
				                          as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMIN surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void PartialCompileTimeSurfaceIMIN(std::size_t      n,
                                       const as3double *B,
				                               as3double       *C,
				                               as3double       *CDerR,
				                               as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMAX surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void PartialCompileTimeSurfaceIMAX(std::size_t      n,
                                       const as3double *B,
				                               as3double       *C,
				                               as3double       *CDerR,
				                               as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMIN surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void PartialCompileTimeSurfaceJMIN(std::size_t      n,
                                       const as3double *B,
				                               as3double       *C,
				                               as3double       *CDerR,
				                               as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMAX surface integration points.
		 *
		 * @param[in] n number of variables processed with the fixed K/M operators.
		 * @param[in] B solution values [n][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [n][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [n][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [n][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void PartialCompileTimeSurfaceJMAX(std::size_t      n,
                                       const as3double *B,
				                               as3double       *C,
				                               as3double       *CDerR,
				                               as3double       *CDerS) override;


		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   volume integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M*M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M*M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M*M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeVolume(const as3double *B,
								           as3double       *C,
								           as3double       *CDerR,
								           as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMIN surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeSurfaceIMIN(const as3double *B,
								                as3double       *C,
								                as3double       *CDerR,
								                as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   IMAX surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeSurfaceIMAX(const as3double *B,
								                as3double       *C,
								                as3double       *CDerR,
								                as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMIN surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeSurfaceJMIN(const as3double *B,
								                as3double       *C,
								                as3double       *CDerR,
								                as3double       *CDerS) override;

		/*!
		 * @brief Function that interpolates the solution and optional parametric derivatives to the
		 *   JMAX surface integration points.
		 *
		 * @param[in] B solution values [N][K*K], with i varying fastest.
		 * @param[out] C optional interpolated values [N][M]; nullptr skips this output.
		 * @param[out] CDerR optional derivative in the r-direction [N][M]; nullptr skips this output.
		 * @param[out] CDerS optional derivative in the s-direction [N][M]; nullptr skips this output.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeSurfaceJMAX(const as3double *B,
								                as3double       *C,
								                as3double       *CDerR,
								                as3double       *CDerS) override;

		/*!
		 * @brief Function that lifts weighted inputs from the volume integration points to the solution
		 *   DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M*M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M*M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M*M]; nullptr omits this
		 *   contribution.
		 * @param[out] C solution residual [N][K*K], overwritten by the volume contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeResidualVolume(const as3double *B,
												           const as3double *BDerR,
												           const as3double *BDerS,
												           as3double       *C) override;

		/*!
		 * @brief Function that lifts weighted inputs from the IMIN surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeResidualSurfaceIMIN(const as3double *B,
												                const as3double *BDerR,
												                const as3double *BDerS,
												                as3double       *C) override;

		/*!
		 * @brief Function that lifts weighted inputs from the IMAX surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeResidualSurfaceIMAX(const as3double *B,
												                const as3double *BDerR,
												                const as3double *BDerS,
												                as3double       *C) override;

		/*!
		 * @brief Function that lifts weighted inputs from the JMIN surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeResidualSurfaceJMIN(const as3double *B,
												                const as3double *BDerR,
												                const as3double *BDerS,
												                as3double       *C) override;

		/*!
		 * @brief Function that lifts weighted inputs from the JMAX surface integration points to the
		 *   solution DOFs.
		 *
		 * @param[in] B optional weighted value contribution [N][M]; nullptr omits this contribution.
		 * @param[in] BDerR optional weighted r-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in] BDerS optional weighted s-derivative contribution [N][M]; nullptr omits this
		 *   contribution.
		 * @param[in,out] C solution residual [N][K*K], incremented by this face contribution.
		 *
		 * @note Arrays are variable-major; the i-index varies fastest within a two-dimensional field.
		 */
		void CompileTimeResidualSurfaceJMAX(const as3double *B,
                                        const as3double *BDerR,
                                        const as3double *BDerS,
                                        as3double       *C) override;

    /*!
     * @brief Function that replaces the element residual by its inverse-mass product.
     *
     * @param[in] minv row-major inverse mass matrix of size [K*K][K*K].
     * @param[in,out] res element residual [N][K*K], replaced by its inverse-mass product.
     * @param[out] tmp scratch storage for at least K*K entries, separate from the residual.
     *
     * @note The scratch buffer must not overlap the residual or inverse mass matrix.
     */
    void CompileTimeApplyInverseMassMatrix(const as3double *minv,
                                           as3double       *res,
                                           as3double       *tmp) override;

	private:
		CMatrixAS3<as3double> mLagrangeInt1D;          ///< 1D Lagrange matrix evaluated at mRInt1D, 
		CMatrixAS3<as3double> mLagrangeInt1DTrans;     ///< Transpose of the 1D Lagrange matrix evaluated at mRInt1D, 
		CMatrixAS3<as3double> mDerLagrangeInt1D;       ///< 1D Derivative of the Lagrange matrix evaluated at mRInt1D,
		CMatrixAS3<as3double> mDerLagrangeInt1DTrans;  ///< Transpose of the 1D derivative of the Lagrange matrix evaluated at mRInt1D,
		CMatrixAS3<as3double> mDerLagrangeMinFace1D;   ///< 1D Derivative of the Lagrange matrix evaluated at: r = -1.0,
		CMatrixAS3<as3double> mDerLagrangeMaxFace1D;   ///< 1D Derivative of the Lagrange matrix evaluated at: r = +1.0,
};

// Definitions of the templated functions.
#include "tensor/tensor_volume.inl"
#include "tensor/tensor_surface_imin.inl"
#include "tensor/tensor_surface_imax.inl"
#include "tensor/tensor_surface_jmin.inl"
#include "tensor/tensor_surface_jmax.inl"
#include "tensor/tensor_residual_volume.inl"
#include "tensor/tensor_residual_surface_imin.inl"
#include "tensor/tensor_residual_surface_imax.inl"
#include "tensor/tensor_residual_surface_jmin.inl"
#include "tensor/tensor_residual_surface_jmax.inl"
#include "tensor/tensor_apply_minv.inl"
