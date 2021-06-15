#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube3d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>

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

      // 2d, Cube, first order
      template< class D >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::cube(2).id(), D, D, 2, 1 >
        : public BDM1Cube2DLocalFiniteElement< D, D >
      {
        static const int numOrientations = 16;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BDM1Cube2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
      };

      // 3d, Cube, first order
      template< class D >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::cube(3).id(), D, D, 3, 1 >
        : public BDM1Cube3DLocalFiniteElement< D, D >
      {
        static const int numOrientations = 64;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BDM1Cube3DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
      };

      // 2d, Cube, second order
      template< class D >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::cube(2).id(), D, D, 2, 2 >
        : public BDM2Cube2DLocalFiniteElement< D, D >
      {
        static const int numOrientations = 16;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BDM2Cube2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
      };


      // 2d, simplex, first order
      template< class D >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::simplex(2).id(), D, D, 2, 1 >
        : public BDM1Simplex2DLocalFiniteElement< D, D >
      {
        static const int numOrientations = 8;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BDM1Simplex2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
      };

      // 2d, simplex, second order
      template< class D >
      struct BDMLocalFiniteElement< Dune::GeometryTypes::simplex(2).id(), D, D, 2, 2 >
        : public BDM2Simplex2DLocalFiniteElement< D, D >
      {
        static const int numOrientations = 8;
        template< class ... Args >
        BDMLocalFiniteElement ( Args && ... args )
          : BDM2Simplex2DLocalFiniteElement< D, D >( std::forward< Args >( args ) ... ) {}
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
