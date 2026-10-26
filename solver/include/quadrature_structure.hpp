///////////////////////////////////////
// TODO: This file needs refactoring...
///////////////////////////////////////
#pragma once

/*!
 * @file
 * @brief Gauss-Jacobi quadrature and its supporting recurrence algorithms.
 */


#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "error_structure.hpp"


/*!
 * @brief Gauss-Jacobi quadrature construction using recurrence and eigensystem routines.
 *
 * @note Classical quadrature helpers follow the IQPACK algorithms of Elhay and Kautsky, implemented in C++ by
 *   John Burkardt.
 */
class CGaussJacobiQuadrature 
{
	public:

  	/*!
  	 * @brief Function that computes Gauss-Jacobi integration points and weights on [a,b].
  	 *
  	 * @param[in] alpha Parameter in the weighting function (b-x)^alpha*(x-a)^beta
  	 * @param[in] beta Parameter in the weighting function (b-x)^alpha*(x-a)^beta
  	 * @param[in] a Lower bound of the integration interval, usually -1.0.
  	 * @param[in] b Upper bound of the integration interval, usually  1.0.
  	 * @param[in,out] GJPoints Location of the Gauss-Jacobi integration points.
  	 * @param[in,out] GJWeights Weights of the Gauss-Jacobi integration points.
  	 *
  	 * @pre GJPoints and GJWeights have the same size, between 1 and 100; their size determines the
  	 *   quadrature order.
  	 */
  	static void GetQuadraturePoints(const double         alpha,    const double         beta,
  	                                const double         a,        const double         b,
  	                                std::vector<double> &GJPoints, std::vector<double> &GJWeights);

	private:
  	// Constructor of the class, disabled.
  	CGaussJacobiQuadrature(void)  = delete;
  	// Destructor of the class, disabled.
  	~CGaussJacobiQuadrature(void) = delete;
		// Disable default copy constructor.
		CGaussJacobiQuadrature(const CGaussJacobiQuadrature&) = delete;
		// Disable default copy operator.
		CGaussJacobiQuadrature& operator=(CGaussJacobiQuadrature&) = delete;

  	/*!
  	 * @brief Function that constructs a classical quadrature rule on its default interval.
  	 *
  	 * @param[in] nt number of knots or recurrence coefficients.
  	 * @param[in] kind identifier of the classical quadrature family.
  	 * @param[in] alpha first weight-function exponent.
  	 * @param[in] beta second weight-function exponent.
  	 * @param[out] t quadrature knots.
  	 * @param[out] wts quadrature weights.
  	 */
  	static void cdgqf(int nt, int kind, double alpha, double beta, double t[], double wts[]);

  	/*!
  	 * @brief Function that constructs and scales a classical quadrature rule to the requested interval.
  	 *
  	 * @param[in] nt number of knots or recurrence coefficients.
  	 * @param[in] kind identifier of the classical quadrature family.
  	 * @param[in] alpha first weight-function exponent.
  	 * @param[in] beta second weight-function exponent.
  	 * @param[in] a lower integration bound.
  	 * @param[in] b upper integration bound.
  	 * @param[out] t quadrature knots.
  	 * @param[out] wts quadrature weights.
  	 */
  	static void cgqf(int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[]);

  	/*!
  	 * @brief Function that builds Jacobi recurrence coefficients and returns the zeroth moment.
  	 *
  	 * @param[in] kind identifier of the classical quadrature family.
  	 * @param[in] m number of knots or recurrence coefficients.
  	 * @param[in] alpha first weight-function exponent.
  	 * @param[in] beta second weight-function exponent.
  	 * @param[out] aj diagonal coefficients of the Jacobi recurrence matrix.
  	 * @param[out] bj off-diagonal coefficients of the Jacobi recurrence matrix.
  	 *
  	 * @return zeroth moment of the selected quadrature weight function.
  	 */
  	static double class_matrix(int kind, int m, double alpha, double beta, double aj[], double bj[]);

  	/*!
  	 * @brief Function that diagonalizes a symmetric tridiagonal matrix and transforms a vector.
  	 *
  	 * @param[in] n number of knots or recurrence coefficients.
  	 * @param[in,out] d diagonal entries of the tridiagonal matrix.
  	 * @param[in,out] e off-diagonal entries of the tridiagonal matrix.
  	 * @param[in,out] z vector transformed by the eigensystem calculation.
  	 */
  	static void imtqlx(int n, double d[], double e[], double z[]);

  	/*!
  	 * @brief Function that validates the parameters of a classical quadrature rule.
  	 *
  	 * @param[in] kind identifier of the classical quadrature family.
  	 * @param[in] m order of the highest moment whose existence is checked.
  	 * @param[in] alpha first weight-function exponent.
  	 * @param[in] beta second weight-function exponent.
  	 */
  	static void parchk(int kind, int m, double alpha, double beta);

  	/*!
  	 * @brief Function that returns the double-precision machine epsilon used by the quadrature routines.
  	 *
  	 * @return double-precision machine epsilon.
  	 */
  	static double r8_epsilon();

  	/*!
  	 * @brief Function that returns the sign of a double-precision value.
  	 *
  	 * @param[in] x x-coordinate data or scalar argument.
  	 *
  	 * @return -1.0 for negative x, otherwise +1.0.
  	 */
  	static double r8_sign(double x);

  	/*!
  	 * @brief Function that scales quadrature knots and weights to a different interval.
  	 *
  	 * @param[in] nt number of knots or recurrence coefficients.
  	 * @param[in] t quadrature knots.
  	 * @param[in] mlt multiplicity of each knot.
  	 * @param[in] wts quadrature weights.
  	 * @param[in] nwts number of weight entries.
  	 * @param[in] ndx starting index of the weights for each knot.
  	 * @param[out] swts scaled quadrature weights.
  	 * @param[out] st scaled quadrature knots.
  	 * @param[in] kind identifier of the classical quadrature family.
  	 * @param[in] alpha first weight-function exponent.
  	 * @param[in] beta second weight-function exponent.
  	 * @param[in] a lower integration bound.
  	 * @param[in] b upper integration bound.
  	 */
  	static void scqf(int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[], 
				             double st[], int kind, double alpha, double beta, double a, double b);

  	/*!
  	 * @brief Function that constructs quadrature knots and weights from Jacobi recurrence coefficients.
  	 *
  	 * @param[in] nt number of knots or recurrence coefficients.
  	 * @param[in] aj diagonal coefficients of the Jacobi recurrence matrix.
  	 * @param[in,out] bj off-diagonal coefficients of the Jacobi recurrence matrix.
  	 * @param[in] zemu zeroth moment of the quadrature weight function.
  	 * @param[out] t quadrature knots.
  	 * @param[out] wts quadrature weights.
  	 */
  	static void sgqf(int nt, double aj[], double bj[], double zemu, double t[], double wts[]);
};
