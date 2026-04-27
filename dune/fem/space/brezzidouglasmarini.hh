#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH

#if HAVE_DUNE_LOCALFUNCTIONS

#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1cube3d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini1simplex2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2cube2d.hh>
#include <dune/localfunctions/brezzidouglasmarini/brezzidouglasmarini2simplex2d.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

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
        static const int numOrientations = 0;
        //static_assert( AlwaysFalse< DomainField >::value, "BDMLocalFiniteElement not implemented for your choice." );
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

    template< class GridPart, class FunctionSpace, int polOrder = -1 >
    class BrezziDouglasMariniLocalFiniteElementMap
    {
      typedef BrezziDouglasMariniLocalFiniteElementMap< GridPart, FunctionSpace, polOrder > ThisType;
      static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                     "GridPart has more than one geometry type." );

      static const unsigned int topologyId = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId;

    public:
      typedef unsigned int KeyType;

      typedef GridPart GridPartType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      typedef PiolaTransformation< typename EntityType::Geometry, FunctionSpace::dimRange > TransformationType;

      static const int dimLocal = GridPart::dimension;

    protected:
      using DomainType = typename FunctionSpace::DomainType;
      using RangeType  = typename FunctionSpace::RangeType;
      typedef LocalBasisTraits<DomainFieldType, FunctionSpace::dimDomain, DomainType,
                               RangeFieldType,  FunctionSpace::dimRange,  RangeType,
                               typename FunctionSpace::JacobianRangeType > LBT;

      // for dynamic space the pOrder is passed as negative value
      static const bool dynamicOrder = polOrder < 0 ;
    public:
      // type of LocalFinitElement
      typedef typename std::conditional< dynamicOrder,
              // virtual LFE type
              LocalFiniteElementVirtualInterface< LBT >,
              // real implementation type
              Impl::BDMLocalFiniteElement< topologyId, DomainFieldType,
                  RangeFieldType, dimLocal, polOrder > > :: type     LocalFiniteElementType;


      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType   LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType  LocalInterpolationType;

      BrezziDouglasMariniLocalFiniteElementMap ( const GridPart &gridPart, const unsigned int ord )
        : orientation_( gridPart ), order_( dynamicOrder ? ord : polOrder )
      {
        if constexpr ( dynamicOrder )
        {
          if ( order_ == 1 )
            this->template createLFE< 1 >();
          else if ( order_ == 2 )
            this->template createLFE< 2 >();
          else
            DUNE_THROW(NotImplemented,"RaviartThomasLocalFiniteElement not implemented for your choice." );
        }
        else
        {
          for ( auto i : range( size() ) )
            map_[ i ].reset( new LocalFiniteElementType( i ) );
        }
      }

      // this is the same for each order, so take order 1 which exists for all combinations
      static constexpr std::size_t numOrientations = Impl::BDMLocalFiniteElement< topologyId, DomainFieldType,
                  RangeFieldType, dimLocal, 1 >::numOrientations;

      static constexpr std::size_t size () { return numOrientations; }

      int order () const { return order_; }

      template< class Entity >
      int order ( const Entity &entity ) const { return order(); }

      template< class Entity >
      auto operator() ( const Entity &e ) const
      {
        unsigned int orient = orientation_( e );
        return std::tuple< std::size_t, const LocalBasisType&, const LocalInterpolationType& >
        { static_cast< std::size_t >( orient ),
          map( orient ).localBasis(),
          map( orient ).localInterpolation() };
      }

      bool hasCoefficients ( const GeometryType &t ) const
      {
        Dune::GeometryType type( GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId, GridPart::dimension );
        return (type == t);
      }

      const LocalCoefficientsType &localCoefficients ( const GeometryType &type ) const
      {
        return map( 0 ).localCoefficients();
      }

      const GridPartType &gridPart () const { return orientation_.gridPart(); }

    private:
      template <int p>
      void createLFE()
      {
        using LFEImpl = Impl::BDMLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, p >;
        if constexpr ( LFEImpl::numOrientations > 0 )
        {
          using LFEObject = LocalFiniteElementVirtualImp< LFEImpl >;
          for ( auto i : range( size() ) )
          {
            LFEImpl imp( i );
            map_[ i ].reset( new LFEObject( imp ) );
          }
        }
        else
        {
          DUNE_THROW(NotImplemented,"BDMLocalFiniteElement not implemented for your choice." );
        }
      }

      const LocalFiniteElementType& map( const std::size_t i ) const
      {
        assert( map_[ i ] );
        return *map_[ i ];
      }

      UniqueFacetOrientation< GridPartType > orientation_;
      std::array< std::unique_ptr< LocalFiniteElementType >, numOrientations > map_;

      const int order_;
    };


    // BrezziDouglasMariniSpace
    // ------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    using BrezziDouglasMariniSpace
            = LocalFiniteElementSpace< BrezziDouglasMariniLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >,
                                       FunctionSpace, Storage >;

    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    using DynamicBrezziDouglasMariniSpace
            = LocalFiniteElementSpace< BrezziDouglasMariniLocalFiniteElementMap< GridPart, FunctionSpace >,
                                       FunctionSpace, Storage >;

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_HH
