/*!
 * @file
 * @brief Collocated compile-time tensor kernels for apply minv.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of the apply minv templated functions in CCollocatedTensorProduct.
//-----------------------------------------------------------------------------------


/*!
 * @brief Function that replaces the element residual by its inverse-mass product.
 *
 * @tparam K number of solution points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] minv packed inverse-mass diagonal [K*K].
 * @param[in,out] res element residual [N][K*K], replaced by its inverse-mass product.
 * @param[out] tmp unused scratch pointer; nullptr is permitted.
 *
 * @note The inverse mass must use collocated quadrature, including physical metrics.
 */
template<std::size_t K, std::size_t N>
void CCollocatedTensorProduct<K, N>::CompileTimeApplyInverseMassMatrix
(
 const as3double *minv,
 as3double       *res,
 as3double       *tmp
)
 /*
	* Function that scales the residual by the diagonal collocated inverse mass.
	*/
{
	// Total number of volume solution points.
	constexpr std::size_t nDOFs = K*K;

	// Scratch storage is not needed for the diagonal inverse mass.
	(void) tmp;

	for(std::size_t iVar=0; iVar<N; iVar++)
	{
		as3double *output = res + iVar * nDOFs;

#pragma omp simd
		for(std::size_t iRow=0; iRow<nDOFs; iRow++)
		{
			output[iRow] *= minv[iRow];
		}
	}
}
