#ifndef DUNE_FEM_SPACE_LAGRANGESPACE_DOFMAPPERCODE_HH
#define DUNE_FEM_SPACE_LAGRANGESPACE_DOFMAPPERCODE_HH

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

    template< class LagrangePointSetContainer >
    struct LagrangeDofMapperCodeFactory
    {
      explicit LagrangeDofMapperCodeFactory ( const LagrangePointSetContainer &lagrangePointSets )
      : lagrangePointSets_( lagrangePointSets )
      {}

      template< class Field, int dim >
      Fem::DofMapperCode operator() ( const Dune::ReferenceElement< Field, dim > &refElement ) const
      {
        const GeometryType type = refElement.type();
        if( lagrangePointSets_.exists( type ) )
          return Fem::compile( refElement, lagrangePointSets_[ type ] );
        else
          return Fem::DofMapperCode();
      }

    private:
      const LagrangePointSetContainer &lagrangePointSets_;
    };

  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGESPACE_DOFMAPPERCODE_HH
