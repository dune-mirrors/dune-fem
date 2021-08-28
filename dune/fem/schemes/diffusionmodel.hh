#ifndef ELLIPTIC_MODEL_HH
#define ELLIPTIC_MODEL_HH
#warning "Deprecated header, #include <dune/fem/schemes/conservationlawmodel.hh>
#include <dune/fem/schemes/conservationlawmodel.hh>

namespace Dune {
namespace Fem {

  template< class GridPart, int dimDomain, int dimRange=dimDomain, class RangeField = double >
  using DiffusionModel = ConservationLawModel< GridPart, dimDomain, dimRange, RangeField >;

  template< class ModelImpl >
  using DiffusionModelWrapper = ConservationLawModelWrapper< ModelImpl >;

  template< class GridPart, int dimDomain, int dimRange=dimDomain, class RangeField = double >
  using DGDiffusionModel = DGConservationLawModel< GridPart, dimDomain, dimRange, RangeField >;

  template < class ModelImpl >
  using DGDiffusionModelWrapper = DGConservationLawModelWrapper< ModelImpl >;

}}

#endif
