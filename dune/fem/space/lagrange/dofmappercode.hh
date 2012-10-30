#ifndef DUNE_FEM_SPACE_LAGRANGE_DOFMAPPER_HH
#define DUNE_FEM_SPACE_LAGRANGE_DOFMAPPER_HH

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/typeindex.hh>

// dune-fem includes
#include <dune/fem/space/dofmapper/code.hh>
#include <dune/fem/space/dofmapper/compile.hh>


namespace Dune
{

  namespace Fem
  {

    // LagrangeDofMapperCodeFactory
    // ----------------------------

    template< class LagrangePointSetContainer >
    struct LagrangeDofMapperCodeFactory
    {
      explicit LagrangeDofMapperCodeFactory ( const LagrangePointSetContainer &lagrangePointSets )
      : lagrangePointSets_( lagrangePointSets )
      {}

      template< class Field, int dim >
      DofMapperCode operator() ( const ReferenceElement< Field, dim > &refElement ) const
      {
        const GeometryType type = refElement.type();
        if( lagrangePointSets_.exists( type ) )
          return compile( refElement, lagrangePointSets_[ type ] );
        else
          return DofMapperCode();
      }

    private:
      const LagrangePointSetContainer &lagrangePointSets_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_DOFMAPPER_HH
