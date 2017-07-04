#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFEMAP_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFEMAP_HH

#include <array>

// dune-geometry types
#include <dune/geometry/type.hh>

// dune-localfunctions includes
#include <dune/fem/space/brezzidouglasmarini/localfiniteelement.hh>
#include <dune/fem/space/basisfunctionset/piolatransformation.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart, class FunctionSpace, int polOrder >
    class BDMLocalFiniteElementMap
    {
      typedef BDMLocalFiniteElementMap< GridPart, FunctionSpace, polOrder > ThisType;
      static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                     "GridPart has more than one geometry type." );

      static const unsigned int topologyId = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId;
    public:
      typedef std::tuple< > KeyType;

      typedef GridPart GridPartType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef PiolaTransformation< typename EntityType::Geometry, FunctionSpace::dimRange > TransformationType;

      static const int dimLocal = GridPart::dimension;

      typedef BDMLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, polOrder > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      template< class ... Args >
      BDMLocalFiniteElementMap ( const GridPart &gridPart, Args ... args )
      : gridPart_( gridPart )
      {
        for( std::size_t i = 0; i < LocalFiniteElementType::numOrientations; ++i )
          map_[ i ] = LocalFiniteElementType( i );
      }

      static std::size_t size () { return LocalFiniteElementType::numOrientations; }

      int order () const { return polOrder; }

      template< class Entity >
      int order ( const Entity &entity ) const { return order(); }

      template< class Entity >
      std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & > operator() ( const Entity &e ) const
      {
        unsigned char orient = orientation( e );
        return std::make_tuple(
          static_cast< std::size_t >( orient ),
          map_[ orient ].localBasis(),
          map_[ orient ].localInterpolation() );
      }

      bool hasCoefficients ( const GeometryType &t ) const
      {
        Dune::GeometryType type( GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId, GridPart::dimension );
        return (type == t);
      }

      const LocalCoefficientsType& localCoefficients ( const GeometryType &type ) const
      {
        return map_[ 0 ].localCoefficients();
      }

      const GridPartType &gridPart () const { return gridPart_; }

    protected:

      // NOTE: this might be cached in future versions
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
      const GridPartType &gridPart_;
      std::array< LocalFiniteElementType, LocalFiniteElementType::numOrientations > map_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_LOCALFEMAP_HH
