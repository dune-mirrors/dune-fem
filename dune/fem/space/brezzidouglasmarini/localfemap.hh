#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFEMAP_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFEMAP_HH

#include <array>

// dune-geometry types
#include <dune/geometry/type.hh>

// dune-localfunctions includes
#include <dune/localfunctions/rannacherturek.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart, class FunctionSpace, int order >
    class BDMLocalFiniteElementMap
    {
      BDMLocalFiniteElementMap< GridPart, FunctionSpace, order > ThisType;
      static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                     "GridPart has more than one geometry type." );

    public:
      typedef std::tuple< > KeyType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      static const int dimLocal = GridPart::dimension;

      typedef BDMLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, order > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      template< class ... Args >
      BDMLocalFiniteElementMap ( const GridPart &gridPart, Args ... args )
      {
        for( std::size_t i = 0; i < LocalFiniteElementType::numOrientations; ++i )
          map_[ i ] = LocalFiniteElementType( i );
      }

      static std::size_t size () const { return LocalFiniteElementType::numOrientations; }

      int order () const { return order; }

      template< class Entity >
      std::tuple< std::size_t, LocalBasisType, LocalInterpolationType > operator() ( const Entity &e ) const
      {
        unsigned char orient = orientation( entity );
        return std::make_tuple(
          static_cast< std::size_t >( orient ),
          map_[ orient ].localBasis(),
          map_[ orient ].localInterpolation() );
      }

      bool hasCoefficient ( const GeometryType &type ) const
      {
        return type.id() == Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId;
      }

      LocalCoefficientsType coefficient ( const GeometryType &type ) const
      {
        return map_[ 0 ].localCoefficient();
      }

    protected:

      template< class Entity >
      unsigned char orientation ( const Entity &entity ) const
      {
        unsigned char ret = 0;
        auto &idxSet = gridPart().indexSet();
        for( auto intersection : intersections( gridPart(), entity ) )
          if( intersection.neighbor() &&
              ( idxSet.index( entity ) < idxSet.index( intersection.outside() ) ) )
            ret |= 1 << intersection.indexInInside();
        return ret;
      }

    private:
      std::array< LocalFiniteElementType, LocalFiniteElementType::numOrientations > map_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFEMAP_HH
