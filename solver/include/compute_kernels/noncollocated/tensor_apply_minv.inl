/*!
 * @file
 * @brief Compile-time tensor kernels for apply minv.
 */

#include <cstddef>



//-----------------------------------------------------------------------------------
// Implementation of: Res = M^{-1} * Res templated function in CTensorProduct.
//-----------------------------------------------------------------------------------

/*!
 * @brief Function that replaces the element residual by its inverse-mass product.
 *
 * @tparam K number of solution points in 1D.
 * @tparam M number of integration points in 1D.
 * @tparam N number of variables.
 *
 * @param[in] minv row-major inverse mass matrix of size [K*K][K*K].
 * @param[in,out] res element residual [N][K*K], replaced by its inverse-mass product.
 * @param[out] tmp scratch storage for at least K*K entries, separate from the residual.
 *
 * @note The scratch buffer must not overlap the residual or inverse mass matrix.
 */
template<std::size_t K, std::size_t M, std::size_t N>
void CTensorProduct<K, M, N>::CompileTimeApplyInverseMassMatrix
(
 const as3double *minv,
 as3double       *res,
 as3double       *tmp
)
 /*
  *
  */
{
  constexpr std::size_t nDOFs = K*K;

  for(std::size_t iVar=0; iVar<N; iVar++)
  {
    const as3double *input  = res + iVar * nDOFs;
    as3double       *output = res + iVar * nDOFs;

    for(std::size_t iRow=0; iRow<nDOFs; iRow++)
    {
      const as3double *matrix_row = minv + iRow * nDOFs;
      as3double             value = C_ZERO;
     
#pragma omp simd reduction(+:value) 
      for(std::size_t iCol=0; iCol<nDOFs; iCol++)
      {
        value += matrix_row[iCol] * input[iCol];
      }
      tmp[iRow] = value;
    }
    for(std::size_t i=0; i<nDOFs; i++) output[i] = tmp[i];
  }
}
