#pragma once

#include "option_structure.hpp"
#include "standard_element_structure.hpp"



/*!
 * @brief An interface class used for the templated tensor specilizations.
 */
class ITensorProduct
{
	public:

		/*!
		 * @brief Constructor of ITensorProduct, which serves as an interface for the tensor product.
		 */
		ITensorProduct(void) = default;
		
		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		virtual ~ITensorProduct(void) = default;


		/*!
		 * @brief Function that implements a generic (run-time) tensor-product volume implementation.
		 *
		 * @param[in] K dimension of DOFs of input points in 1D.
		 * @param[in] N column dimension of both input and output vector.
		 * @param[in] M dimension of DOFs of output points in 1D.
 		 * @param[in] A pointer to the (transposed) interpolation polynomial in 1D, dimension: [K][M].
 		 * @param[in] ADer pointer to the (transposed) derivative of the interpolation polynomial in 1D, dimension: [K][M].
 		 * @param[in] B pointer to input 2D data, dimension: [K*K][N].
 		 * @param[out] C pointer to output interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerR pointer to output x-derivative of interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerS pointer to output y-derivative of interpolated data, dimension: [M*M][N].
		 */
		void CustomVolume(size_t           K,
				              size_t           N,
 				              size_t           M,
 				              const as3double *A,
 				              const as3double *ADer,
 				              const as3double *B,
 				              as3double       *C,
 				              as3double       *CDerR,
 				              as3double       *CDerS);

		/*!
		 * @brief Function that implements a generic (run-time) tensor-product IMIN surface implementation.
		 *
		 * @param[in] K dimension of DOFs of input points in 1D.
		 * @param[in] N column dimension of both input and output vector.
		 * @param[in] M dimension of DOFs of output points in 1D.
 		 * @param[in] A pointer to the (transposed) interpolation polynomial in 1D, dimension: [K][M].
 		 * @param[in] ADer pointer to the (transposed) derivative of the interpolation polynomial in 1D, dimension: [K][M].
		 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial on the IMIN surface, dimension: [K].
 		 * @param[in] B pointer to input 2D data, dimension: [K*K][N].
 		 * @param[out] C pointer to output interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerR pointer to output x-derivative of interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerS pointer to output y-derivative of interpolated data, dimension: [M*M][N].
		 */
		void CustomSurfaceIMIN(size_t           K,
				                   size_t           N,
 				                   size_t           M,
 				                   const as3double *A,
 				                   const as3double *ADer,
													 const as3double *ADerFace,
 				                   const as3double *B,
 				                   as3double       *C,
 				                   as3double       *CDerR,
 				                   as3double       *CDerS);

		/*!
		 * @brief Function that implements a generic (run-time) tensor-product IMAX surface implementation.
		 *
		 * @param[in] K dimension of DOFs of input points in 1D.
		 * @param[in] N column dimension of both input and output vector.
		 * @param[in] M dimension of DOFs of output points in 1D.
 		 * @param[in] A pointer to the (transposed) interpolation polynomial in 1D, dimension: [K][M].
 		 * @param[in] ADer pointer to the (transposed) derivative of the interpolation polynomial in 1D, dimension: [K][M].
		 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial on the IMAX surface, dimension: [K].
 		 * @param[in] B pointer to input 2D data, dimension: [K*K][N].
 		 * @param[out] C pointer to output interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerR pointer to output x-derivative of interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerS pointer to output y-derivative of interpolated data, dimension: [M*M][N].
		 */
		void CustomSurfaceIMAX(size_t           K,
				                   size_t           N,
 				                   size_t           M,
 				                   const as3double *A,
 				                   const as3double *ADer,
													 const as3double *ADerFace,
 				                   const as3double *B,
 				                   as3double       *C,
 				                   as3double       *CDerR,
 				                   as3double       *CDerS);

		/*!
		 * @brief Function that implements a generic (run-time) tensor-product JMIN surface implementation.
		 *
		 * @param[in] K dimension of DOFs of input points in 1D.
		 * @param[in] N column dimension of both input and output vector.
		 * @param[in] M dimension of DOFs of output points in 1D.
 		 * @param[in] A pointer to the (transposed) interpolation polynomial in 1D, dimension: [K][M].
 		 * @param[in] ADer pointer to the (transposed) derivative of the interpolation polynomial in 1D, dimension: [K][M].
		 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial on the JMIN surface, dimension: [K].
 		 * @param[in] B pointer to input 2D data, dimension: [K*K][N].
 		 * @param[out] C pointer to output interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerR pointer to output x-derivative of interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerS pointer to output y-derivative of interpolated data, dimension: [M*M][N].
		 */
		void CustomSurfaceJMIN(size_t           K,
				                   size_t           N,
 				                   size_t           M,
 				                   const as3double *A,
 				                   const as3double *ADer,
													 const as3double *ADerFace,
 				                   const as3double *B,
 				                   as3double       *C,
 				                   as3double       *CDerR,
 				                   as3double       *CDerS);

		/*!
		 * @brief Function that implements a generic (run-time) tensor-product JMAX surface implementation.
		 *
		 * @param[in] K dimension of DOFs of input points in 1D.
		 * @param[in] N column dimension of both input and output vector.
		 * @param[in] M dimension of DOFs of output points in 1D.
 		 * @param[in] A pointer to the (transposed) interpolation polynomial in 1D, dimension: [K][M].
 		 * @param[in] ADer pointer to the (transposed) derivative of the interpolation polynomial in 1D, dimension: [K][M].
		 * @param[in] ADerFace pointer to the derivative of the interpolation polynomial on the JMAX surface, dimension: [K].
 		 * @param[in] B pointer to input 2D data, dimension: [K*K][N].
 		 * @param[out] C pointer to output interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerR pointer to output x-derivative of interpolated data, dimension: [M*M][N].
 		 * @param[out] CDerS pointer to output y-derivative of interpolated data, dimension: [M*M][N].
		 */
		void CustomSurfaceJMAX(size_t           K,
				                   size_t           N,
 				                   size_t           M,
 				                   const as3double *A,
 				                   const as3double *ADer,
													 const as3double *ADerFace,
 				                   const as3double *B,
 				                   as3double       *C,
 				                   as3double       *CDerR,
 				                   as3double       *CDerS);


		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product on the volume integration points.
		 */
		virtual void Volume(const size_t     N,
				                const as3double *B,
								        as3double       *C,
								        as3double       *CDerR,
								        as3double       *CDerS) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product on the IMIN surface integration points.
		 */
		virtual void SurfaceIMIN(const size_t     N,
				                     const as3double *B,
								             as3double       *C,
								             as3double       *CDerR,
								             as3double       *CDerS) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product on the IMAX surface integration points.
		 */
		virtual void SurfaceIMAX(const size_t     N,
				                     const as3double *B,
								             as3double       *C,
								             as3double       *CDerR,
								             as3double       *CDerS) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product on the JMIN surface integration points.
		 */
		virtual void SurfaceJMIN(const size_t     N,
				                     const as3double *B,
								             as3double       *C,
								             as3double       *CDerR,
								             as3double       *CDerS) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product on the JMAX surface integration points.
		 */
		virtual void SurfaceJMAX(const size_t     N,
				                     const as3double *B,
								             as3double       *C,
								             as3double       *CDerR,
								             as3double       *CDerS) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product that computes the residual from the volume terms. 
		 */
		virtual void ResidualVolume(const size_t     N,
				                        const as3double *B,
																const as3double *BDerR,
																const as3double *BDerS,
																as3double       *C) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product that computes the residual from the IMIN surface terms. 
		 */
		virtual void ResidualSurfaceIMIN(const size_t     N,
				                             const as3double *B,
																     const as3double *BDerR,
																     const as3double *BDerS,
																     as3double       *C) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product that computes the residual from the IMAX surface terms. 
		 */
		virtual void ResidualSurfaceIMAX(const size_t     N,
				                             const as3double *B,
																     const as3double *BDerR,
																     const as3double *BDerS,
																     as3double       *C) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product that computes the residual from the JMIN surface terms. 
		 */
		virtual void ResidualSurfaceJMIN(const size_t     N,
				                             const as3double *B,
																     const as3double *BDerR,
																     const as3double *BDerS,
																     as3double       *C) = 0;

		/*!
		 * @brief Pure virtual function for a (compile-time) tensor product that computes the residual from the JMAX surface terms. 
		 */
		virtual void ResidualSurfaceJMAX(const size_t     N,
				                             const as3double *B,
																     const as3double *BDerR,
																     const as3double *BDerS,
																     as3double       *C) = 0;

	protected:

	private:
		// Disable default copy constructor.
		ITensorProduct(const ITensorProduct&) = delete;
		// Disable default copy operator.
		ITensorProduct& operator=(ITensorProduct&) = delete;	
};

//-----------------------------------------------------------------------------------

/*!
 * @brief A templated implementation for different tensor-product functions, based on a fixed (K,M).
 */
template<size_t K, size_t M>
class CTensorProduct final: public ITensorProduct
{
	public:
		
		/*!
		 * @brief Constructor which initializes this specialized tensor-product class.
		 *
		 * @param[in] standard_element pointer to the standard element of this specialization.
		 */
		CTensorProduct(CStandardElement *standard_element) :
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

		// Default destructor.
		~CTensorProduct(void) final = default;


		/*!
		 * @brief Function that computes the tensor product on the volume integration points.
		 *
		 * @param[in] N number of variables.
		 * @param[in] B pointer to the input values at the solution DOFs.
		 * @param[out] C pointer to the interpolated values at the integration points.
		 * @param[out] CDerR pointer to the derivative in the r-direction at the integration points.
		 * @param[out] CDerS pointer to the derivative in the s-direction at the integration points.
		 */
		void Volume(const size_t     N,
				        const as3double *B,
								as3double       *C,
								as3double       *CDerR,
								as3double       *CDerS) final;

		/*!
		 * @brief Function that computes the tensor product on the surface integration points of IMIN.
		 *
		 * @param[in] N number of variables.
		 * @param[in] B pointer to the input values at the solution DOFs.
		 * @param[out] C pointer to the interpolated values at the integration points.
		 * @param[out] CDerR pointer to the derivative in the r-direction at the integration points.
		 * @param[out] CDerS pointer to the derivative in the s-direction at the integration points.
		 */
		void SurfaceIMIN(const size_t     N,
				             const as3double *B,
								     as3double       *C,
								     as3double       *CDerR,
								     as3double       *CDerS) final;

		/*!
		 * @brief Function that computes the tensor product on the surface integration points of IMAX.
		 *
		 * @param[in] N number of variables.
		 * @param[in] B pointer to the input values at the solution DOFs.
		 * @param[out] C pointer to the interpolated values at the integration points.
		 * @param[out] CDerR pointer to the derivative in the r-direction at the integration points.
		 * @param[out] CDerS pointer to the derivative in the s-direction at the integration points.
		 */
		void SurfaceIMAX(const size_t     N,
				             const as3double *B,
								     as3double       *C,
								     as3double       *CDerR,
								     as3double       *CDerS) final;

		/*!
		 * @brief Function that computes the tensor product on the surface integration points of JMIN.
		 *
		 * @param[in] N number of variables.
		 * @param[in] B pointer to the input values at the solution DOFs.
		 * @param[out] C pointer to the interpolated values at the integration points.
		 * @param[out] CDerR pointer to the derivative in the r-direction at the integration points.
		 * @param[out] CDerS pointer to the derivative in the s-direction at the integration points.
		 */
		void SurfaceJMIN(const size_t     N,
				             const as3double *B,
								     as3double       *C,
								     as3double       *CDerR,
								     as3double       *CDerS) final;

		/*!
		 * @brief Function that computes the tensor product on the surface integration points of JMAX.
		 *
		 * @param[in] N number of variables.
		 * @param[in] B pointer to the input values at the solution DOFs.
		 * @param[out] C pointer to the interpolated values at the integration points.
		 * @param[out] CDerR pointer to the derivative in the r-direction at the integration points.
		 * @param[out] CDerS pointer to the derivative in the s-direction at the integration points.
		 */
		void SurfaceJMAX(const size_t     N,
				             const as3double *B,
								     as3double       *C,
								     as3double       *CDerR,
								     as3double       *CDerS) final;

		/*!
		 * @brief Function for a (compile-time) tensor product that computes the residual from the volume terms. 
		 *        Note, this also resets the residual.
		 */
		void ResidualVolume(const size_t     N,
		                    const as3double *B,
												const as3double *BDerR,
												const as3double *BDerS,
												as3double       *C) final;

		/*!
		 * @brief Function for a (compile-time) tensor product that computes the residual from the IMIN surface terms. 
		 */
		void ResidualSurfaceIMIN(const size_t     N,
		                         const as3double *B,
												     const as3double *BDerR,
												     const as3double *BDerS,
												     as3double       *C) final;

		/*!
		 * @brief Function for a (compile-time) tensor product that computes the residual from the IMAX surface terms. 
		 */
		void ResidualSurfaceIMAX(const size_t     N,
		                         const as3double *B,
												     const as3double *BDerR,
												     const as3double *BDerS,
												     as3double       *C) final;

		/*!
		 * @brief Function for a (compile-time) tensor product that computes the residual from the JMIN surface terms. 
		 */
		void ResidualSurfaceJMIN(const size_t     N,
		                         const as3double *B,
												     const as3double *BDerR,
												     const as3double *BDerS,
												     as3double       *C) final;

		/*!
		 * @brief Function for a (compile-time) tensor product that computes the residual from the JMAX surface terms. 
		 */
		void ResidualSurfaceJMAX(const size_t     N,
		                         const as3double *B,
												     const as3double *BDerR,
												     const as3double *BDerS,
												     as3double       *C) final;

	private:
		CMatrixAS3<as3double> mLagrangeInt1D;         ///< 1D Lagrange matrix evaluated at mRInt1D, 
		CMatrixAS3<as3double> mLagrangeInt1DTrans;    ///< Transpose of the 1D Lagrange matrix evaluated at mRInt1D, 
		CMatrixAS3<as3double> mDerLagrangeInt1D;      ///< 1D Derivative of the Lagrange matrix evaluated at mRInt1D,
		CMatrixAS3<as3double> mDerLagrangeInt1DTrans; ///< Transpose of the 1D derivative of the Lagrange matrix evaluated at mRInt1D,
		CMatrixAS3<as3double> mDerLagrangeMinFace1D;  ///< 1D Derivative of the Lagrange matrix evaluated at: r = -1.0,
		CMatrixAS3<as3double> mDerLagrangeMaxFace1D;  ///< 1D Derivative of the Lagrange matrix evaluated at: r = +1.0,
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

