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

  template< class LagrangePointSetContainer >
  struct LagrangeDofMapperCodeFactory
  {
    explicit LagrangeDofMapperCodeFactory ( const LagrangePointSetContainer &lagrangePointSets )
    : lagrangePointSets_( lagrangePointSets )
    {}

    template< class Field, int dim >
    Fem::DofMapperCode operator() ( const GenericReferenceElement< Field, dim > &refElement ) const
    {
      std::size_t gtIdx = LocalGeometryTypeIndex::index( refElement.type() );
      if( lagrangePointSets_[ gtIdx ] )
        return Fem::compile( refElement, *lagrangePointSets_[ gtIdx ] );
      else
        return Fem::DofMapperCode();
    }

  private:
    const LagrangePointSetContainer &lagrangePointSets_;
  };

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGESPACE_DOFMAPPERCODE_HH
