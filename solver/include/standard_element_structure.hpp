#pragma once 

/*!
 * @file
 * @brief Reference-element nodes, quadrature and one-dimensional basis operators.
 */


#include <cstddef>

#include <cmath>
#include "option_structure.hpp"
#include "config_structure.hpp"
#include "quadrature_structure.hpp"
#include "polynomial_structure.hpp"
#include "integration_rule_structure.hpp"


/*!
 * @brief Reference-element nodes, quadrature and basis operators shared by physical elements of the same
 *   discretization.
 */
class CStandardElement 
{
	public:

		/*!
		 * @brief Constructor that initializes reference nodes, quadrature weights and polynomial
		 *   operators.
		 *
		 * @param[in] typeSolutionDistribution distribution of the solution nodes.
		 * @param[in] nPolySolution polynomial degree of the solution in 1D.
		 * @param[in] nBonusQuadratureFactor non-negative nb in the nodal-DG target degree (2 + nb)*p.
		 * @param[in] typeQuadrature quadrature rule on the reference interval.
		 * @param[in] typeDGMethod DG formulation; DG-SEM requires zero bonus and LGL nodes/quadrature.
		 */
		CStandardElement(ETypeDOF        typeSolutionDistribution,
				             std::size_t     nPolySolution,
				             std::size_t     nBonusQuadratureFactor,
                     ETypeQuadrature typeQuadrature,
                     ETypeDGMethod   typeDGMethod);

		/*!
		 * @brief Destructor that releases the resources owned by this object.
		 */
		~CStandardElement(void);


		/*!
		 * @brief Getter function which returns the solution polynomial order.
		 *
		 * @return mNPolySol.
		 */
		std::size_t GetnPolySol(void) const {return mNPolySol;}

		/*!
		 * @brief Getter function which returns the number of integration points in one reference
		 *   direction.
		 *
		 * @return number of 1D points supplied by the integration rule.
		 */
		std::size_t GetnInt1D(void) const {return mIntegrationRuleContainer->GetnInt1D();}

		/*!
		 * @brief Getter function which returns the number of solution points in one reference
		 *   direction.
		 *
		 * @return mNSol1D.
		 */
		std::size_t GetnSol1D(void) const {return mNSol1D;}

		/*!
		 * @brief Getter function which returns the number of integration points in 2D.
		 *
		 * @return number of tensor-product points supplied by the integration rule.
		 */
		std::size_t GetnInt2D(void) const {return mIntegrationRuleContainer->GetnInt2D();}

		/*!
		 * @brief Getter function which returns the number of solution points in 2D.
		 *
		 * @return mNSol2D.
		 */
		std::size_t GetnSol2D(void) const {return mNSol2D;}

		/*!
		 * @brief Getter function which returns the Lagrange interpolation matrix in 1D, at the
		 *   integration points.
		 *
		 * @return mLagrangeInt1D.
		 */
		const CMatrixAS3<as3double> &GetLagrangeInt1D(void) const {return mLagrangeInt1D;}

		/*!
		 * @brief Getter function which returns the differentiation matrix in 1D, at the integration
		 *   points.
		 *
		 * @return mDerLagrangeInt1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeInt1D(void) const {return mDerLagrangeInt1D;}

		/*!
		 * @brief Getter function which returns the differentiation matrix in 1D, at the solution
		 *   points.
		 *
		 * @return mDerLagrangeSol1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeSol1D(void) const {return mDerLagrangeSol1D;}

		/*!
		 * @brief Getter function which returns the transposed differentiation matrix in 1D, at the
		 *   solution points.
		 *
		 * @return mDerLagrangeSol1DTrans.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeSol1DTrans(void) const {return mDerLagrangeSol1DTrans;}

		/*!
		 * @brief Getter function which returns the transposed Lagrange interpolation matrix in 1D, at
		 *   the integration points.
		 *
		 * @return mLagrangeInt1DTrans.
		 */
		const CMatrixAS3<as3double> &GetLagrangeInt1DTrans(void) const {return mLagrangeInt1DTrans;}

		/*!
		 * @brief Getter function which returns the transposed differentiation matrix in 1D, at the
		 *   integration points.
		 *
		 * @return mDerLagrangeInt1DTrans.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeInt1DTrans(void) const {return mDerLagrangeInt1DTrans;}

		/*!
		 * @brief Getter function which returns the derivative of the Lagrange polynomial in 1D, at the
		 *   MIN surface.
		 *
		 * @return mDerLagrangeMinFace1D.
		 */
		const CMatrixAS3<as3double> &GetDerLagrangeMinFace1D(void) const {return mDerLagrangeMinFace1D;}

		/*!
		 * @brief Getter function which returns the derivative of the Lagrange polynomial in 1D, at the
		 *   MAX surface.
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
		ETypeDOF GetTypeDOFsSol(void) const {return mTypeDOFsSol;}

		/*!
		 * @brief Getter function which returns the integration points in 1D, over the reference
		 *   element.
		 *
		 * @return integration-rule nodes on [-1,+1].
		 */
		const CMatrixAS3<as3double> &GetrInt1D(void) const {return mIntegrationRuleContainer->GetrInt1D();}

		/*!
		 * @brief Getter function which returns the integration weights in 1D, over the reference
		 *   element.
		 *
		 * @return integration-rule weights on [-1,+1].
		 */
		const CMatrixAS3<as3double> &GetwInt1D(void) const {return mIntegrationRuleContainer->GetwInt1D();}

		/*!
		 * @brief Getter function which returns the integration weights in 2D, over the reference
		 *   element.
		 *
		 * @return tensor-product integration weights on [-1,+1]^2.
		 */
		const CMatrixAS3<as3double> &GetwInt2D(void) const {return mIntegrationRuleContainer->GetwInt2D();}

		/*!
		 * @brief Getter function which returns the shared reference-element integration rule.
		 *
		 * @return non-owning pointer to the rule; valid for the lifetime of this element.
		 */
		const CIntegrationRule* GetIntegrationRule(void) const {return mIntegrationRuleContainer.get();}

	protected:

	private:
		ETypeDOF      mTypeDOFsSol;  ///< Type of solution DOFs.
		std::size_t   mNPolySol;     ///< Solution polynomial order.
		std::size_t   mNSol1D;       ///< Number of solution nodes in 1D.
		std::size_t   mNSol2D;       ///< Number of solution nodes in 2D.

		std::unique_ptr<CIntegrationRule> mIntegrationRuleContainer; ///< Owned quadrature nodes and 1D/2D weights.

		CMatrixAS3<as3double> mRSol1D;  ///< Solution points on the reference element in 1D.

		CMatrixAS3<as3double> mVandermondeSol1D;      ///< 1D Vandermonde matrix evaluated at mRSol1D:,
																						      ///< [mNSol1D][mNSol1D].
		CMatrixAS3<as3double> mVandermondeInt1D;      ///< 1D Vandermonde matrix evaluated at GetrInt1D(),
																						      ///< [GetnInt1D()][mNSol1D].

		CMatrixAS3<as3double> mDerVandermondeSol1D;   ///< 1D Derivative of the Vandermonde matrix evaluated at mRSol1D, 
																						      ///< with dimensions: [mNSol1D][mNSol1D].
		CMatrixAS3<as3double> mDerVandermondeInt1D;   ///< 1D Derivative of the Vandermonde matrix evaluated at GetrInt1D(),
																						      ///< with dimensions: [GetnInt1D()][mNSol1D].

		CMatrixAS3<as3double> mLagrangeInt1D;         ///< 1D Lagrange matrix evaluated at GetrInt1D(),
																								  ///< with dimensions: [GetnInt1D()][mNSol1D].
		CMatrixAS3<as3double> mLagrangeInt1DTrans;    ///< Transpose of the 1D Lagrange matrix evaluated at GetrInt1D(),
																								  ///< with dimensions: [mNSol1D][GetnInt1D()].
    
		CMatrixAS3<as3double> mDerLagrangeSol1D;      ///< 1D Derivative of the Lagrange matrix evaluated at mRSol1D,
																								  ///< with dimensions: [mNSol1D][mNSol1D].
		CMatrixAS3<as3double> mDerLagrangeSol1DTrans; ///< Transpose of the 1D derivative of the Lagrange matrix evaluated at mRSol1D,
																								  ///< with dimensions: [mNSol1D][mNSol1D].

		CMatrixAS3<as3double> mDerLagrangeInt1D;      ///< 1D Derivative of the Lagrange matrix evaluated at GetrInt1D(),
																								  ///< with dimensions: [GetnInt1D()][mNSol1D].
		CMatrixAS3<as3double> mDerLagrangeInt1DTrans; ///< Transpose of the 1D derivative of the Lagrange matrix evaluated at GetrInt1D(),
																								  ///< with dimensions: [mNSol1D][GetnInt1D()].

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
                               std::size_t           nDOFs1D,
                               CMatrixAS3<as3double> &rDOFs1D);

		// Disable default constructor.
		CStandardElement(void) = delete;
		// Disable default copy constructor.
		CStandardElement(const CStandardElement&) = delete;
		// Disable default copy operator.
		CStandardElement& operator=(CStandardElement&) = delete;	
};
