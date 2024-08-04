#pragma once

#include <random>
#include "option_structure.hpp"
#include "gtest/gtest.h"


namespace NTest_LA
{
	// Function that initializes a matrix sequentially.
	void InitMatrixSeq(CMatrixAS3<as3double> &A);

	// Function that initializes a matrix randomly.
	void InitMatrixRand(CMatrixAS3<as3double> &A);

	// Function that computes the correct matrix transpose.
	void MatrixTranspose(CMatrixAS3<as3double> &A, 
			                 CMatrixAS3<as3double> &B);

	// Function that computes the correct matrix-matrix multiplication.
	void MatrixMult(CMatrixAS3<as3double> &A, 
			            CMatrixAS3<as3double> &B, 
									CMatrixAS3<as3double> &C);

	// Function that checks the L1-norm of the matrix differences.
	void MatrixErrorNormLinf(CMatrixAS3<as3double> &A, 
			                     CMatrixAS3<as3double> &B,
												   as3double              tol);

	// Function that computes the inverse of a 3x3 matrix.
	void InverseMatrix3x3(CMatrixAS3<as3double> &A,
			                  CMatrixAS3<as3double> &B);

	// Helper function for debugging, by printing a matrix.
	void PrintAS3Matrix(CMatrixAS3<as3double> &A, 
			                std::string            g,
											unsigned short         p = 6);
}
