///////////////////////////////////////
// TODO: This file needs refactoring...
///////////////////////////////////////
#pragma once

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "error_log.hpp"


// Definition of the class CGaussJacobiQuadrature.
class CGaussJacobiQuadrature 
{
	public:

  	/*!
  	 * \brief Function, which serves as the API to compute the integration points
  	          and weights.
  	 * \param[in]     alpha     Parameter in the weighting function (b-x)^alpha*(x-a)^beta
  	                            in the Gauss Jacobi rule.
  	 * \param[in]     beta      Parameter in the weighting function (b-x)^alpha*(x-a)^beta
  	                            in the Gauss Jacobi rule.
  	 * \param[in]     a         Lower bound of the integration interval, usually -1.0.
  	 * \param[in]     b         Upper bound of the integration interval, usually  1.0.
  	 * \param[in,out] GJPoints  Location of the Gauss-Jacobi integration points.
  	 * \param[in,out] GJWeights Weights of the Gauss-Jacobi integration points.
  	 */
  	static void GetQuadraturePoints(const double         alpha,    const double         beta,
  	                                const double         a,        const double         b,
  	                                std::vector<double> &GJPoints, std::vector<double> &GJWeights);

	private:
  	// Constructor of the class, disabled.
  	CGaussJacobiQuadrature()  = delete;
  	// Destructor of the class, disabled.
  	~CGaussJacobiQuadrature() = delete;

		/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static void cdgqf(int nt, int kind, double alpha, double beta, double t[], double wts[]);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static void cgqf(int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[]);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static double class_matrix(int kind, int m, double alpha, double beta, double aj[], double bj[]);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static void imtqlx(int n, double d[], double e[], double z[]);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static void parchk(int kind, int m, double alpha, double beta);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static double r8_epsilon();

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static double r8_sign(double x);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static void scqf(int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[], 
				             double st[], int kind, double alpha, double beta, double a, double b);

  	/*!
  	 * \brief Function in the original implementation of John Burkardt to compute
  	          the integration points of the Gauss-Jacobi quadrature rule.
  	 */
  	static void sgqf(int nt, double aj[], double bj[], double zemu, double t[], double wts[]);
};
