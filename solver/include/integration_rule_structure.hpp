#pragma once

/*!
 * @file
 * @brief Reference quadrature nodes, 1D/2D weights and one-dimensional polynomial exactness.
 */


#include <cstddef>

#include "option_structure.hpp"
#include "data_structure.hpp"


/*!
 * @brief Quadrature rule on [-1,+1], with tensor-product weights on [-1,+1]^2.
 *
 * @note Points and weights are generated together and cannot be modified
 *   through the public interface.
 */
class CIntegrationRule
{
	public:

		/*!
		 * @brief Constructor that initializes the quadrature points and weights.
		 *
		 * @param[in] typeQuadrature quadrature rule.
		 * @param[in] nPolyIntegrand minimum polynomial degree to integrate exactly in 1D.
		 */
		CIntegrationRule(ETypeQuadrature typeQuadrature,
		                 std::size_t     nPolyIntegrand);

		/*!
		 * @brief Default destructor.
		 */
		~CIntegrationRule(void) = default;


		/*!
		 * @brief Function that computes the point count required for a polynomial degree.
		 *
		 * @param[in] typeQuadrature quadrature rule on the reference interval.
		 * @param[in] nPolyIntegrand minimum polynomial degree to integrate exactly in 1D.
		 *
		 * @return required number of integration points in 1D, without allocating quadrature data.
		 */
		static std::size_t ComputePointCount1D(ETypeQuadrature typeQuadrature,
		                                       std::size_t     nPolyIntegrand);

		/*!
		 * @brief Getter function which returns the quadrature rule.
		 *
		 * @return mTypeQuadrature.
		 */
		ETypeQuadrature GetTypeQuadrature(void) const
		{
			return mTypeQuadrature;
		}

		/*!
		 * @brief Getter function which returns the number of integration
		 *   points in 1D.
		 *
		 * @return number of entries in mWInt1D.
		 */
		std::size_t GetnInt1D(void) const
		{
			return mWInt1D.size();
		}

		/*!
		 * @brief Getter function which returns the number of integration
		 *   points in 2D.
		 *
		 * @return number of entries in mWInt2D.
		 */
		std::size_t GetnInt2D(void) const
		{
			return mWInt2D.size();
		}

		/*!
		 * @brief Getter function which returns the integration points in 1D.
		 *
		 * @return mRInt1D.
		 */
		const CMatrixAS3<as3double> &GetrInt1D(void) const
		{
			return mRInt1D;
		}

		/*!
		 * @brief Getter function which returns the integration weights in 1D.
		 *
		 * @return mWInt1D.
		 */
		const CMatrixAS3<as3double> &GetwInt1D(void) const
		{
			return mWInt1D;
		}

		/*!
		 * @brief Getter function which returns the integration weights in 2D.
		 *
		 * @return mWInt2D.
		 */
		const CMatrixAS3<as3double> &GetwInt2D(void) const
		{
			return mWInt2D;
		}

		/*!
		 * @brief Getter function which returns the highest polynomial degree
		 *   integrated exactly by this quadrature rule.
		 *
		 * @note Exactness refers to polynomial integration with unit weight
		 *   on [-1,+1], in exact arithmetic.
		 */
		std::size_t GetnPolyExact1D(void) const;

		/*!
		 * @brief Function that checks whether all polynomials up to the
		 *   specified degree are integrated exactly.
		 *
		 * @param[in] nPoly polynomial degree of the integrand.
		 *
		 * @return true if nPoly does not exceed the degree of exactness.
		 */
		bool IsExactForPolynomial1D(std::size_t nPoly) const
		{
			return nPoly <= GetnPolyExact1D();
		}

    /*!
     * @brief Equality operator that checks whether two objects represent
     *   the same quadrature rule.
     *
     * @param[in] other integration rule to compare against.
     *
     * @return true if the quadrature types and point counts are identical.
     */
    bool operator==(const CIntegrationRule &other) const
    {
    	return mTypeQuadrature == other.mTypeQuadrature && GetnInt1D() == other.GetnInt1D();
    }


    /*!
     * @brief Inequality operator that checks whether two objects represent
     *   different quadrature rules.
     *
     * @param[in] other integration rule to compare against.
     *
     * @return true if the quadrature rules differ.
     */
    bool operator!=(const CIntegrationRule &other) const
    {
    	return !(*this == other);
    }

	private:

		ETypeQuadrature       mTypeQuadrature; ///< Quadrature rule.
		CMatrixAS3<as3double> mRInt1D;         ///< Integration points in 1D.
		CMatrixAS3<as3double> mWInt1D;         ///< Integration weights in 1D.
    CMatrixAS3<as3double> mWInt2D;         ///< Integration weights in 2D.
		
    /*!
		 * @brief Function that initializes Gauss-Legendre points and weights.
		 *
		 * @param[in] nPolyIntegrand minimum polynomial degree to integrate exactly in 1D.
		 */
		void InitializeGaussLegendre(std::size_t nPolyIntegrand);

		/*!
		 * @brief Function that initializes Legendre-Gauss-Lobatto points
		 *   and weights.
		 *
		 * @param[in] nPolyIntegrand minimum polynomial degree to integrate exactly in 1D.
		 */
		void InitializeLegendreGaussLobatto(std::size_t nPolyIntegrand);


		// Disable default constructor.
		CIntegrationRule(void) = delete;
		// Disable default copy constructor.
		CIntegrationRule(const CIntegrationRule&) = delete;
		// Disable default copy operator.
		CIntegrationRule& operator=(const CIntegrationRule&) = delete;
};
