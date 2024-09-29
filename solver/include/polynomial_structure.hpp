#pragma once 

#include "option_structure.hpp"
#include "linear_algebra.hpp"

/*!
 * @brief A namespace used for storing polynomial functionalities.
 */
namespace NPolynomialUtility
{
  /*!
	 * @brief Function that computes the integration rule in 1D.
	 *
	 * @param[in] npoly input polynomial order.
	 *
	 * @return number of integration nodes in 1D.
	 */
  unsigned short IntegrationRule1D(unsigned short npoly);

  /*!
	 * @brief Function that computes the Lagrange basis functions in 1D via 2 Vandermonde matrices: L = V1*inv(V2).
	 *
	 * @param[in] V1 Vandermonde matrix in 1D at evaluation points.
	 * @param[in] V2 Vandermonde matrix in 1D at solution points.
	 * @param[out] L1D Lagrange interpolation matrix in 1D.
	 */
	void LagrangeBasis1D(CMatrixAS3<as3double> &V1,
											 CMatrixAS3<as3double> &V2,
											 CMatrixAS3<as3double> &L1D);

  /*!
	 * @brief Function that computes the Lagrange differentiation matrix in 1D via 2 Vandermonde matrices: dL = dV1*inv(V2).
	 *
	 * @param[in] V1 derivative of the Vandermonde matrix in 1D at evaluation points.
	 * @param[in] V2 Vandermonde matrix in 1D at solution points.
	 * @param[out] dL1D Lagrange differentiation matrix in 1D.
	 */
	void DerivativeLagrangeBasis1D(CMatrixAS3<as3double> &dV1,
											           CMatrixAS3<as3double> &V2,
											           CMatrixAS3<as3double> &dL1D);

  /*!
	 * @brief Function that computes the derivative of the Lagrange basis function at the min/max faces in 1D.
	 *
	 * @param[in] dLSol1D differentiation (square) matrix in 1D, evaluated at the solution points.
	 * @param[out] dLmin derivative of the Lagrange function at the min face: r = -1.
	 * @param[out] dLmax derivative of the Lagrange function at the max face: r = +1.
	 */
	void DerivativeLagrangeBasisFace1D(const CMatrixAS3<as3double> &dLSol1D,
											               CMatrixAS3<as3double>       &dLmin,
											               CMatrixAS3<as3double>       &dLmax);

  /*!
	 * @brief Function that computes the the Legendre basis functions in 1D.
	 *
	 * @param[in] nCol number of columns corresponding to different orders.
	 * @param[in] rDOFs1D vector of DOFs in 1D.
	 * @param[out] rBasis1D vector containing the Legendre basis in 1D.
	 */
	void OrthonormalLegendreBasis1D(const unsigned int           nCol, 
			                            const CMatrixAS3<as3double> &rDOFs1D,
																	CMatrixAS3<as3double>       &rBasis1D);

  /*!
	 * @brief Function that computes the derivative of the Legendre basis functions in 1D.
	 *
	 * @param[in] nCol number of columns corresponding to different orders.
	 * @param[in] rDOFs1D vector of DOFs in 1D.
	 * @param[out] rDerBasis1D vector containing the derivative of the Legendre basis in 1D.
	 */
	void DerivativeOrthonormalLegendreBasis1D(const unsigned int           nCol, 
			                                      const CMatrixAS3<as3double> &rDOFs1D,
																	          CMatrixAS3<as3double>       &rDerBasis1D);

  /*!
	 * @brief Function that computes the normalized Legendre polynomial coefficient at a point.
	 *
	 * @param[in] n order of the polynomial.
	 * @param[in] r evaluation point over the interval: [-1,+1].
	 * 
	 * @return value of the Legendre polynomial at the input point and order.
	 */
	as3double NormalizedLegendrePolynomial(const int n, const as3double r);

  /*!
	 * @brief Function that computes the derivative of the normalized Legendre polynomial coefficient at a point.
	 *
	 * @param[in] n order of the polynomial.
	 * @param[in] r evaluation point over the interval: [-1,+1].
	 * 
	 * @return value of the derivative of the Legendre polynomial at the input point and order.
	 */
	as3double DerivativeNormalizedLegendrePolynomial(const int n, const as3double r);



  /*!
	 * @brief Function that computes the value of the normalized Jacobi polynomial at a point.
	 *
	 * @param[in] n order of the function.
	 * @param[in] alpha polynomial coefficient.
	 * @param[in] beta polynomial coefficient.
	 * @param[in] r evaluation point over the interval: [-1,+1].
	 * 
	 * @return value of the normalized Jacobi polynomial based on input order, point and parameters.
	 */
	as3double NormalizedJacobiPolynomial(const int       n, 
			                                 const int       alpha,
																		   const int       beta,
																		   const as3double r);


  /*!
	 * @brief Function that computes the value of the derivative of the normalized Jacobi polynomial at a point.
	 *
	 * @param[in] n order of the function.
	 * @param[in] alpha polynomial coefficient.
	 * @param[in] beta polynomial coefficient.
	 * @param[in] r evaluation point over the interval: [-1,+1].
	 * 
	 * @return value of the derivative of the normalized Jacobi polynomial based on input order, point and parameters.
	 */
	as3double DerivativeNormalizedJacobiPolynomial(const int       n, 
			                                           const int       alpha,
																		             const int       beta,
																		             const as3double r);


}
