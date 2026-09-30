

//-----------------------------------------------------------------------------------
// Implementation of: Res = M^{-1} * Res templated function in CTensorProduct.
//-----------------------------------------------------------------------------------

template<size_t K, size_t M, size_t N>
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
  constexpr size_t nDOFs = K*K;

  for(size_t iVar=0; iVar<N; iVar++)
  {
    const as3double *input  = res + iVar * nDOFs;
    as3double       *output = res + iVar * nDOFs;

    for(size_t iRow=0; iRow<nDOFs; iRow++)
    {
      const as3double *matrix_row = minv + iRow * nDOFs;
      as3double             value = C_ZERO;
     
#pragma omp simd reduction(+:value) 
      for(size_t iCol=0; iCol<nDOFs; iCol++)
      {
        value += matrix_row[iCol] * input[iCol];
      }
      tmp[iRow] = value;
    }
    for(size_t i=0; i<nDOFs; i++) output[i] = tmp[i];
  }
}
