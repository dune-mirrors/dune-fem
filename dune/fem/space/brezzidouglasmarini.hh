#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/localfunctions/brezzidouglasmarini.hh>

#include <dune/fem/space/common/uniquefacetorientation.hh>
#include <dune/fem/space/basisfunctionset/piolatransformation.hh>
#include <dune/fem/space/localfiniteelement/space.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

      // BDMLocalFiniteElement
      // ---------------------

      template< unsigned int id, class DomainField, class RangeField, int dimension, int order >
      struct BDMLocalFiniteElement
      {
        static_assert( AlwaysFalse< DomainField >::value, "BDMLocalFiniteElement not implemented for your choice." );
      };

      // The following local finite elements are implemented

      // cube
      template< class D, int dimension, int order >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::cube(dimension).id(), D, D, dimension, order >
        : public BrezziDouglasMariniCubeLocalFiniteElement< D, D, dimension, order >
      {
        // 2d order 1 and 2, 3d order 1
        static_assert( order >= 1 && order <= 4-dimension, "Selected order not available");
        typedef BrezziDouglasMariniCubeLocalFiniteElement< D, D, dimension, order > BaseType;
        static const int numOrientations = dimension == 2 ? 16 : 64;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BaseType( std::forward< Args >( args ) ... ) {}
      };

      // simplex
      template< class D, int dimension, int order  >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::simplex(dimension).id(), D, D, dimension, order >
        : public BrezziDouglasMariniSimplexLocalFiniteElement< D, D, dimension, order >
      {
        // 2d order 1 and 2, 3d order 1
        static_assert( order >= 1 && order <= 4-dimension, "Selected order not available");

        typedef BrezziDouglasMariniSimplexLocalFiniteElement< D, D, dimension, order > BaseType;
        static const int numOrientations = dimension == 2 ? 8 : 16;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BaseType( std::forward< Args >( args ) ... ) {}
      };
    }

    // BrezziDouglasMariniLocalFiniteElementMap
    // ----------------------------------------

    template< class GridPart, class FunctionSpace, int polOrder >
    class BrezziDouglasMariniLocalFiniteElementMap
    {
      typedef BrezziDouglasMariniLocalFiniteElementMap< GridPart, FunctionSpace, polOrder > ThisType;
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

      typedef Impl::BDMLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, polOrder > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType   LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType  LocalInterpolationType;

      template< class ... Args >
      BrezziDouglasMariniLocalFiniteElementMap ( const GridPart &gridPart, Args ... args )
        : orientation_( gridPart )
      {
        for( std::size_t i = 0; i < LocalFiniteElementType::numOrientations; ++i )
          map_[ i ] = LocalFiniteElementType( i );
      }

      static std::size_t size () { return LocalFiniteElementType::numOrientations; }

      int order () const { return polOrder; }

      template< class Entity >
      int order ( const Entity &entity ) const { return order(); }

      template< class Entity >
      auto operator() ( const Entity &e ) const
      {
        unsigned int orient = orientation_( e );
        return std::tuple< std::size_t, const LocalBasisType&, const LocalInterpolationType& >
        { static_cast< std::size_t >( orient ),
          map_[ orient ].localBasis(),
          map_[ orient ].localInterpolation() };
      }

      bool hasCoefficients ( const GeometryType &t ) const
      {
        Dune::GeometryType type( GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId, GridPart::dimension );
        return (type == t);
      }

      const LocalCoefficientsType &localCoefficients ( const GeometryType &type ) const
      {
        return map_[ 0 ].localCoefficients();
      }

      const GridPartType &gridPart () const { return orientation_.gridPart(); }

    private:
      UniqueFacetOrientation< GridPartType > orientation_;
      std::array< LocalFiniteElementType, LocalFiniteElementType::numOrientations > map_;
    };


    // BrezziDouglasMariniSpace
    // ------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    using BrezziDouglasMariniSpace
            = LocalFiniteElementSpace< BrezziDouglasMariniLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >,
                                       FunctionSpace, Storage >;

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
