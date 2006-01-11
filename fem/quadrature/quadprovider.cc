#ifndef QUADRATURE_QUADPROVIDER_CC
#define QUADRATURE_QUADPROVIDER_CC

#include "quadprovider.hh"

namespace Dune {
  template <typename ct>
  std::vector<CubeQuadrature<ct, 1>*> QuadratureProvider<ct, 1>::
  quads_(CubeQuadrature<ct, 1>::maxOrder(), 0);

  template <typename ct>
  std::vector<SimplexQuadrature<ct, 2>*> QuadratureProvider<ct, 2>::
  triangleQuads_(SimplexQuadrature<ct, 2>::maxOrder(), 0);

  template <typename ct>
  std::vector<CubeQuadrature<ct, 2>*> QuadratureProvider<ct, 2>::
  quadrilateralQuads_(CubeQuadrature<ct, 2>::maxOrder(), 0);

  template <typename ct>
  std::vector<SimplexQuadrature<ct, 3>*> QuadratureProvider<ct, 3>::
  tetraQuads_(SimplexQuadrature<ct, 3>::maxOrder(), 0);
  
  template <typename ct>
  std::vector<CubeQuadrature<ct, 3>*> QuadratureProvider<ct, 3>::
  hexaQuads_(CubeQuadrature<ct, 3>::maxOrder(), 0);
  
  template <typename ct>
  std::vector<PrismQuadrature<ct>*> QuadratureProvider<ct, 3>::
  prismQuads_(PrismQuadrature<ct>::maxOrder(), 0);
  
  template <typename ct>
  std::vector<PyramidQuadrature<ct>*> QuadratureProvider<ct, 3>::
  pyramidQuads_(PyramidQuadrature<ct>::maxOrder(), 0);
}

#endif
