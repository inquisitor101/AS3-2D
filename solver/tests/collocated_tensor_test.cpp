/*!
 * @file
 * @brief Compare collocated kernels with general tensors and polynomial/adjoint identities.
 */
#include "compute_kernels/TensorProducts.hpp"
#include "data/CMatrixAS3.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace
{
using Vec = std::vector<as3double>;
void Check(const Vec &a, const Vec &b)
{
  if(a.size() != b.size()) throw std::runtime_error("Size mismatch.");
  for(std::size_t i=0; i<a.size(); i++)
    if(!std::isfinite(a[i]) || !std::isfinite(b[i]) ||
       std::abs(a[i]-b[i]) > 2.e-11*(1.0+std::abs(b[i])))
      throw std::runtime_error("Collocated kernel mismatch at entry " + std::to_string(i));
}

using Eval = void (ITensorProduct::*)(const as3double*, as3double*, as3double*, as3double*);
using Partial = void (ITensorProduct::*)(std::size_t, const as3double*, as3double*, as3double*, as3double*);
using Lift = void (ITensorProduct::*)(const as3double*, const as3double*, const as3double*, as3double*);
constexpr Eval eval[] = {&ITensorProduct::CompileTimeVolume,
  &ITensorProduct::CompileTimeSurfaceIMIN, &ITensorProduct::CompileTimeSurfaceIMAX,
  &ITensorProduct::CompileTimeSurfaceJMIN, &ITensorProduct::CompileTimeSurfaceJMAX};
constexpr Partial partial[] = {&ITensorProduct::PartialCompileTimeVolume,
  &ITensorProduct::PartialCompileTimeSurfaceIMIN, &ITensorProduct::PartialCompileTimeSurfaceIMAX,
  &ITensorProduct::PartialCompileTimeSurfaceJMIN, &ITensorProduct::PartialCompileTimeSurfaceJMAX};
constexpr Lift lift[] = {&ITensorProduct::CompileTimeResidualVolume,
  &ITensorProduct::CompileTimeResidualSurfaceIMIN, &ITensorProduct::CompileTimeResidualSurfaceIMAX,
  &ITensorProduct::CompileTimeResidualSurfaceJMIN, &ITensorProduct::CompileTimeResidualSurfaceJMAX};

template<std::size_t K, std::size_t N>
void Run()
{
  CStandardElement element(ETypeNodalDistribution::LGL, K-1, 0,
      ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, ETypeDGMethod::DG_SEM);
  CCollocatedTensorProduct<K,N> collocated(&element);
  CTensorProduct<K,K,N> general(&element);
  ITensorProduct &a=collocated, &b=general;
  if(a.GetnSol1D()!=K || a.GetnInt1D()!=K || a.GetnVar()!=N)
    throw std::runtime_error("Incorrect dimensions.");
  // Partial calls must support n independent of N, including n > N.
  const std::size_t n=N+2;
  Vec field(n*K*K);
  for(std::size_t i=0; i<field.size(); i++) field[i]=std::sin(0.37*static_cast<double>(i))+0.2;

  for(std::size_t op=0; op<5; op++)
  {
    const std::size_t points=op ? K : K*K;
    for(unsigned mask=0; mask<8; mask++)
    {
      Vec av(n*points, -13.), ar=av, as=av, bv=av, br=av, bs=av;
      (a.*eval[op])(field.data(), mask&1 ? av.data():nullptr,
          mask&2 ? ar.data():nullptr, mask&4 ? as.data():nullptr);
      (b.*eval[op])(field.data(), mask&1 ? bv.data():nullptr,
          mask&2 ? br.data():nullptr, mask&4 ? bs.data():nullptr);
      Check(av,bv); Check(ar,br); Check(as,bs);
      (a.*partial[op])(n, field.data(), mask&1 ? av.data():nullptr,
          mask&2 ? ar.data():nullptr, mask&4 ? as.data():nullptr);
      (b.*partial[op])(n, field.data(), mask&1 ? bv.data():nullptr,
          mask&2 ? br.data():nullptr, mask&4 ? bs.data():nullptr);
      Check(av,bv); Check(ar,br); Check(as,bs);
      (a.*partial[op])(0, nullptr, nullptr, nullptr, nullptr);

      Vec f(N*points), r=f, s=f;
      for(std::size_t i=0; i<f.size(); i++)
      {
        f[i]=std::cos(0.19*static_cast<double>(i));
        r[i]=std::sin(0.41*static_cast<double>(i));
        s[i]=std::cos(0.73*static_cast<double>(i));
      }
      Vec ac(N*K*K, 3.25), bc=ac;
      (a.*lift[op])(mask&1 ? f.data():nullptr, mask&2 ? r.data():nullptr,
          mask&4 ? s.data():nullptr, ac.data());
      (b.*lift[op])(mask&1 ? f.data():nullptr, mask&2 ? r.data():nullptr,
          mask&4 ? s.data():nullptr, bc.data());
      Check(ac,bc);
      // Independently check that lifting is the adjoint of evaluation.
      Vec value(N*points), dr=value, ds=value;
      (a.*eval[op])(field.data(), value.data(), dr.data(), ds.data());
      double lhs=0.0, rhs=0.0;
      for(std::size_t i=0; i<f.size(); i++)
        lhs+=(mask&1 ? value[i]*f[i]:0.)+(mask&2 ? dr[i]*r[i]:0.)+(mask&4 ? ds[i]*s[i]:0.);
      for(std::size_t i=0; i<ac.size(); i++) rhs+=field[i]*(ac[i]-(op ? 3.25:0.));
      Check(Vec{lhs},Vec{rhs});
    }
  }

  // Differentiate r*s exactly; face derivatives use reference signs at both ends.
  const auto &x=element.GetrSol1D();
  Vec polynomial(N*K*K), value(N*K*K), dr=value, ds=value, expectedR=value, expectedS=value;
  for(std::size_t v=0; v<N; v++)
    for(std::size_t j=0; j<K; j++)
      for(std::size_t i=0; i<K; i++)
      {
        const auto q=v*K*K+j*K+i;
        polynomial[q]=x[i]*x[j]; expectedR[q]=x[j]; expectedS[q]=x[i];
      }
  a.CompileTimeVolume(polynomial.data(), value.data(), dr.data(), ds.data());
  Check(value,polynomial); Check(dr,expectedR); Check(ds,expectedS);

  Vec diagonal(K*K), minv(K*K*K*K, 0.), residual(N*K*K, 2.), expected=residual, scratch(K*K);
  for(std::size_t q=0; q<K*K; q++)
  {
    diagonal[q]=1./element.GetwInt2D()[q];
    minv[q*(K*K+1)]=diagonal[q];
    for(std::size_t v=0; v<N; v++) expected[v*K*K+q]*=minv[q*(K*K+1)];
  }
  Vec reference(N*K*K, 2.);
  b.CompileTimeApplyInverseMassMatrix(minv.data(), reference.data(), scratch.data());
  a.CompileTimeApplyInverseMassMatrix(diagonal.data(), residual.data(), nullptr);
  Check(residual,reference); Check(residual,expected);
  std::cout << "Passed K=" << K << " N=" << N << '\n';
}
}

int main(int argc, char **argv)
{
  try
  {
    if(argc>1)
    {
      const std::string test = argv[1];
      if(test=="null") { CCollocatedTensorProduct<6,4> tensor(nullptr); }
      if(test=="wrong_size")
      {
        CStandardElement element(ETypeNodalDistribution::LGL, 4, 0,
            ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, ETypeDGMethod::DG_SEM);
        CCollocatedTensorProduct<6,4> tensor(&element);
        return 1;
      }
      if(test=="no_min_endpoint" || test=="no_max_endpoint")
      {
        CStandardElement element(ETypeNodalDistribution::LGL, 5, 0,
            ETypeQuadrature::LEGENDRE_GAUSS_LOBATTO, ETypeDGMethod::DG_SEM);
        // Simulate matching node sets without one endpoint. The constructor must
        // reject this fixture before using its differentiation operator.
        const std::size_t i = test=="no_min_endpoint" ? 0 : 5;
        const as3double coordinate = i==0 ? -0.9 : 0.9;
        const_cast<CMatrixAS3<as3double>&>(element.GetrSol1D())[i] = coordinate;
        const_cast<CMatrixAS3<as3double>&>(element.GetrInt1D())[i] = coordinate;
        CCollocatedTensorProduct<6,4> tensor(&element);
        return 1;
      }
      CStandardElement element(ETypeNodalDistribution::LGL, 5, 0,
          ETypeQuadrature::GAUSS_LEGENDRE, ETypeDGMethod::NODAL_DG);
      CCollocatedTensorProduct<6,4> tensor(&element);
      return 1;
    }
    Run<2,1>(); Run<3,4>(); Run<4,4>(); Run<5,4>(); Run<6,4>();
    Run<7,4>(); Run<8,4>(); Run<9,4>(); Run<10,4>();
  }
  catch(const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
