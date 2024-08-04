#pragma once 

#include <cmath>
#include "option_structure.hpp"
#include "config_structure.hpp"
#include "quadrature_structure.hpp"
#include "polynomial_structure.hpp"


/*!
 * @brief A class used for initializing and definining all operations on a standard element in reference space.
 */
class CStandardElement 
{
	public:
		/*!
		 * @brief Default constructor of CStandardElement, which initializes a CStandardElement class.
		 *
		 * @param[in] config_container configuration/dictionary container.
		 * @param[in] iZone current zone ID.
		 */
		CStandardElement(CConfig        *config_container,
						         unsigned short  iZone);

		/*!
		 * @brief Destructor, which frees any allocated memory.
		 */
		~CStandardElement(void);



		/*!
		 * @brief Getter function which returns the solution polynomial order.
		 *
		 * @return mNPolySol.
		 */
		size_t GetnPolySol(void) const {return static_cast<size_t>(mNPolySol);}

		/*!
		 * @brief Getter function which returns the number of integration points in 1D.
		 *
		 * @return mNInt1D.
		 */
		size_t GetnInt1D(void) const {return static_cast<size_t>(mNInt1D);}

		/*!
		 * @brief Getter function which returns the number of solution points in 1D.
		 *
		 * @return mNSol1D.
		 */
		size_t GetnSol1D(void) const {return static_cast<size_t>(mNSol1D);}

		/*!
		 * @brief Getter function which returns the number of integration points in 2D.
		 *
		 * @return mNInt2D.
		 */
		size_t GetnInt2D(void) const {return static_cast<size_t>(mNInt2D);}

		/*!
		 * @brief Getter function which returns the number of solution points in 2D.
		 *
		 * @return mNSol2D.
		 */
		size_t GetnSol2D(void) const {return static_cast<size_t>(mNSol2D);}

		/*!
		 * @brief Getter function which returns the Lagrange interpolation matrix in 1D, at the integration points.
		 *
		 * @return mLagrangeInt1D.
		 */
		const CMatrixAS3<as3double> &GetLagrangeInt1D(void) const {return mLagrangeInt1D;}

		/*!
		 * @brief Getter function which returns the differentiation matrix in 1D, at the integration points.
		 *
		 * @return mDerLagrangeInt1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeInt1D(void) const {return mDerLagrangeInt1D;}

		/*!
		 * @brief Getter function which returns the differentiation matrix in 1D, at the solution points.
		 *
		 * @return mDerLagrangeSol1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeSol1D(void) const {return mDerLagrangeSol1D;}

		/*!
		 * @brief Getter function which returns the transposed differentiation matrix in 1D, at the solution points.
		 *
		 * @return mDerLagrangeSol1DTrans.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeSol1DTrans(void) const {return mDerLagrangeSol1DTrans;}

		/*!
		 * @brief Getter function which returns the transposed Lagrange interpolation matrix in 1D, at the integration points.
		 *
		 * @return mLagrangeInt1DTrans.
		 */
		const CMatrixAS3<as3double> &GetLagrangeInt1DTrans(void) const {return mLagrangeInt1DTrans;}

		/*!
		 * @brief Getter function which returns the transposed differentiation matrix in 1D, at the integration points.
		 *
		 * @return mDerLagrangeInt1DTrans.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeInt1DTrans(void) const {return mDerLagrangeInt1DTrans;}

		/*!
		 * @brief Getter function which returns the derivative of the Lagrange polynomial in 1D, at the MIN surface.
		 *
		 * @return mDerLagrangeMinFace1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeMinFace1D(void) const {return mDerLagrangeMinFace1D;}

		/*!
		 * @brief Getter function which returns the derivative of the Lagrange polynomial in 1D, at the MAX surface.
		 *
		 * @return mDerLagrangeMaxFace1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeMaxFace1D(void) const {return mDerLagrangeMaxFace1D;}

		/*!
		 * @brief Getter function which returns the solution DOFs in 1D, over the reference element.
		 *
		 * @return mRSol1D.
		 */
		const CMatrixAS3<as3double> &GetrSol1D(void) const {return mRSol1D;}

		/*!
		 * @brief Getter function which returns the type of DOFs.
		 *
		 * @return mTypeDOFsSol.
		 */
		const ETypeDOF GetTypeDOFsSol(void) const {return mTypeDOFsSol;}

		/*!
		 * @brief Getter function which returns the integration points in 1D, over the reference element.
		 *
		 * @return mRInt1D.
		 */
		const CMatrixAS3<as3double> &GetrInt1D(void) const {return mRInt1D;}

		/*!
		 * @brief Getter function which returns the integration weights in 1D, over the reference element.
		 *
		 * @return mWInt1D.
		 */
		const CMatrixAS3<as3double> &GetwInt1D(void) const {return mWInt1D;}

		/*!
		 * @brief Getter function which returns the integration weights in 2D, over the reference element.
		 *
		 * @return mWInt2D.
		 */
		const CMatrixAS3<as3double> &GetwInt2D(void) const {return mWInt2D;}



	protected:

	private:
		unsigned short        mNPolySol;    ///< Solution polynomial order.
		unsigned short        mNPolyGrid;   ///< Grid element polynomial order.
		ETypeDOF              mTypeDOFsSol; ///< Type of solution DOFs.
		unsigned int          mNSol1D;      ///< Number of solution nodes in 1D.
		unsigned int          mNSol2D;      ///< Number of solution nodes in 2D.
		unsigned int          mNInt1D;      ///< Number of integration nodes in 1D.
		unsigned int          mNInt2D;      ///< Number of integration nodes in 2D.
		
		CMatrixAS3<as3double> mRSol1D;      ///< Solution points on the reference element in 1D.
    CMatrixAS3<as3double> mRInt1D;      ///< Integration nodes on the reference element in 1D.
    CMatrixAS3<as3double> mWInt1D;      ///< Integration weights on the reference element in 1D.
		CMatrixAS3<as3double> mWInt2D;      ///< Integration weights on the reference element in 2D.

		CMatrixAS3<as3double> mVandermondeSol1D;      ///< 1D Vandermonde matrix evaluated at mRSol1D:,
																						      ///< [mNSol1D][mNSol1D].
		CMatrixAS3<as3double> mVandermondeInt1D;      ///< 1D Vandermonde matrix evaluated at mRInt1D,
																						      ///< [mNInt1D][mNSol1D].

		CMatrixAS3<as3double> mDerVandermondeSol1D;   ///< 1D Derivative of the Vandermonde matrix evaluated at mRSol1D, 
																						      ///< with dimensions: [mNSol1D][mNSol1D].
		CMatrixAS3<as3double> mDerVandermondeInt1D;   ///< 1D Derivative of the Vandermonde matrix evaluated at mRInt1D,
																						      ///< with dimensions: [mNInt1D][mNSol1D].

		CMatrixAS3<as3double> mLagrangeInt1D;         ///< 1D Lagrange matrix evaluated at mRInt1D, 
																								  ///< with dimensions: [mNInt1D][mNSol1D].
		CMatrixAS3<as3double> mLagrangeInt1DTrans;    ///< Transpose of the 1D Lagrange matrix evaluated at mRInt1D, 
																								  ///< with dimensions: [mNSol1D][mNInt1D].
    
		CMatrixAS3<as3double> mDerLagrangeSol1D;      ///< 1D Derivative of the Lagrange matrix evaluated at mRSol1D,
																								  ///< with dimensions: [mNSol1D][mNSol1D].
		CMatrixAS3<as3double> mDerLagrangeSol1DTrans; ///< Transpose of the 1D derivative of the Lagrange matrix evaluated at mRSol1D,
																								  ///< with dimensions: [mNSol1D][mNSol1D].

		CMatrixAS3<as3double> mDerLagrangeInt1D;      ///< 1D Derivative of the Lagrange matrix evaluated at mRInt1D,
																								  ///< with dimensions: [mNInt1D][mNSol1D].
		CMatrixAS3<as3double> mDerLagrangeInt1DTrans; ///< Transpose of the 1D derivative of the Lagrange matrix evaluated at mRInt1D,
																								  ///< with dimensions: [mNSol1D][mNInt1D].

		CMatrixAS3<as3double> mDerLagrangeMinFace1D;  ///< 1D Derivative of the Lagrange matrix evaluated at: r = -1.0,
																								  ///< with dimensions: [mNSol1D].
		CMatrixAS3<as3double> mDerLagrangeMaxFace1D;  ///< 1D Derivative of the Lagrange matrix evaluated at: r = +1.0,
																								  ///< with dimensions: [mNSol1D].


    /*!
		 * @brief Function that computes the location of the DOFs in 1D.
		 *
		 * @param[in] typeDOFs type of DOFs.
		 * @param[in] nDOFs1D number of DOFs in 1D.
		 * @param[out] rDOFs1D vector of DOFs in 1D over a reference element.
		 */
    void ComputeLocationDOFs1D(ETypeDOF               typeDOFs,
                               unsigned int           nDOFs1D,
                               CMatrixAS3<as3double> &rDOFs1D);

    /*!
		 * @brief Function that computes the integration properties over the reference element.
		 *
		 * @param[in] nInt1D number of integration points in 1D.
		 * @param[out] rInt1D integration points in 1D.
		 * @param[out] wInt1D integration weights in 1D.
		 * @param[out] wInt2D integration weights in 2D.
		 */
    void InitializeQuadrature(unsigned int           nInt1D,
                              CMatrixAS3<as3double> &rInt1D,
                              CMatrixAS3<as3double> &wInt1D,
                              CMatrixAS3<as3double> &wInt2D);


   /*!
		 * @brief Function that computes the integration rule.
		 *
		 * @param[in] npoly input polynomial order.
		 *
		 * @return number of integration nodes in 1D.
		 */
    unsigned int IntegrationRule(unsigned short npoly);



};
