#ifndef DUNE_FEM_QUADRATUREIMP_INLINE_HH
#define DUNE_FEM_QUADRATUREIMP_INLINE_HH

#include "quadratureimp.hh"

namespace Dune
{
  
  template <class ct, int dim>
  inline TestQuadrature<ct, dim>::TestQuadrature(const GeometryType& geo, int order) :
    QuadratureImp<ct, dim>(IdProvider::instance().newId()),
    geo_(geo),
    order_(order)
  {}

  template <class ct, int dim>
  inline void TestQuadrature<ct, dim>::
  newQuadraturePoint(const CoordinateType& c, ct w) 
  {
    this->addQuadraturePoint(c, w);
  }

} // end namespace Dune

#endif
