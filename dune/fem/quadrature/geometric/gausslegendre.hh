#ifndef DUNE_FEM_QUADRATURE_GEOMETRIC_GAUSSLEGENDRE_HH
#define DUNE_FEM_QUADRATURE_GEOMETRIC_GAUSSLEGENDRE_HH

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include "geometry.hh"

namespace Dune
{

  namespace Fem
  {

    // GaussLegendreQuadrature
    // -----------------------

    template< class Field, int dim >
    class GaussLegendreQuadrature
    : public GeometryQuadrature< Dune::QuadratureRule< Field, dim > >
    {
      using BaseType = GeometryQuadrature< Dune::QuadratureRule< Field, dim > >;

    public:
      GaussLegendreQuadrature ( Dune::GeometryType type, int order )
        : BaseType( Dune::QuadratureRules< Field, dim >::rule( type, order, Dune::QuadratureType::GaussLegendre ) )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_QUADRATURE_GEOMETRIC_GAUSSLEGENDRE_HH
