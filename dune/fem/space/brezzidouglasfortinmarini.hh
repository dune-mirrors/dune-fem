#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASFORTINMARINI_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASFORTINMARINI_HH

#if HAVE_DUNE_LOCALFUNCTIONS

// C++ includes
#include <array>
#include <tuple>

// dune-common includes
#include <dune/common/math.hh>
#include <dune/common/typetraits.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-localfunction includes
#include <dune/localfunctions/brezzidouglasfortinmarini/bdfmcube.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

// dune-fem includes
#include <dune/fem/space/common/uniquefacetorientation.hh>
#include <dune/fem/space/basisfunctionset/piolatransformation.hh>
#include <dune/fem/space/localfiniteelement/space.hh>


namespace Dune
{

  namespace Fem
  {

    // BrezziDouglasFortinMariniLocalFiniteElementMap
    // ----------------------------------------------

    template< class GridPart, class FunctionSpace, int polOrder = -1 >
    class BrezziDouglasFortinMariniLocalFiniteElementMap
    {
      using hasSingleGeometryType = GridPartCapabilities::hasSingleGeometryType< GridPart >;

      static_assert( hasSingleGeometryType::v, "`GridPart` has more the one geometry type." );
      static_assert( GridPart::dimension == 2, "`BDFM` basis only defined for quadrilaterals." );

      static constexpr int dimLocal = GridPart::dimension;
      static constexpr unsigned int topologyId = hasSingleGeometryType::topologyId;
      static_assert( topologyId == Dune::GeometryTypes::cube(dimLocal).id(),
                     "Only defined for cube grids." );

      static_assert( dimLocal == FunctionSpace::dimRange, "`dimRange` has to be equal to `GridPart::dimension`" );

      using Geometry = typename GridPart::template Codim< 0 >::EntityType::Geometry;

      // for dynamic space the pOrder is passed as negative value
      static const bool dynamicOrder = polOrder < 0 ;
    public:
      using GridPartType = GridPart;

      using DomainFieldType = typename FunctionSpace::DomainFieldType;
      using RangeFieldType = typename FunctionSpace::RangeFieldType;

    protected:
      using DomainType = typename FunctionSpace::DomainType;
      using RangeType  = typename FunctionSpace::RangeType;
      typedef LocalBasisTraits<DomainFieldType, FunctionSpace::dimDomain, DomainType,
                               RangeFieldType,  FunctionSpace::dimRange,  RangeType,
                               typename FunctionSpace::JacobianRangeType > LBT;

    public:
      // type of LocalFiniteElement
      typedef typename std::conditional< dynamicOrder,
              // virtual LFE type
              LocalFiniteElementVirtualInterface< LBT >,
              // real implementation type
              BDFMCubeLocalFiniteElement< DomainFieldType, RangeFieldType, dimLocal, polOrder > > :: type     LocalFiniteElementType;

      typedef unsigned int KeyType;

      using TransformationType = PiolaTransformation< Geometry, FunctionSpace::dimRange >;

      using LocalBasisType          = typename LocalFiniteElementType::Traits::LocalBasisType;
      using LocalCoefficientsType   = typename LocalFiniteElementType::Traits::LocalCoefficientsType;
      using LocalInterpolationType  = typename LocalFiniteElementType::Traits::LocalInterpolationType;

      static constexpr auto size () -> std::size_t { return Dune::power( int(2), int(2*dimLocal) ); }

      BrezziDouglasFortinMariniLocalFiniteElementMap ( const GridPartType& gridPart, const unsigned int ord )
        : orientation_( gridPart ), order_( dynamicOrder ? ord : polOrder )
      {
        if constexpr ( dynamicOrder )
        {
          if ( order_ == 1 )
            this->template createLFE< 1 >();
          else if ( order_ == 2 )
            this->template createLFE< 2 >();
          else if ( order_ == 3 )
            this->template createLFE< 3 >();
          else
            DUNE_THROW(NotImplemented,"BDFMLocalFiniteElement not implemented for your choice." );
        }
        else
        {
          for ( auto i : range( size() ) )
            map_[ i ].reset( new LocalFiniteElementType( i ) );
        }
      }

      int order () const { return order_; }

      template< class Entity >
      int order ( const Entity& entity ) const { return order(); }

      template< class Entity >
      auto operator() ( const Entity& entity ) const
        -> std::tuple< std::size_t, const LocalBasisType&, const LocalInterpolationType& >
      {
        auto o = orientation_( entity );
        return std::tuple< std::size_t, const LocalBasisType&, const LocalInterpolationType& >(
            static_cast< std::size_t >( o ), map( o ).localBasis(), map( o ).localInterpolation() );
      }

      auto localCoefficients ( const GeometryType& type ) const
        -> const LocalCoefficientsType&
      {
        return map( 0 ).localCoefficients();
      }

      bool hasCoefficients ( const GeometryType& type ) const { return type == GeometryType( topologyId, dimLocal ); }
      auto gridPart () const -> const GridPartType& { return orientation_.gridPart(); }

    private:
      template <int p>
      void createLFE()
      {
        using LFEImpl = BDFMCubeLocalFiniteElement< DomainFieldType, RangeFieldType, dimLocal, p >;
        //if constexpr ( LFEImpl::numOrientations > 0 )
        {
          using LFEObject = LocalFiniteElementVirtualImp< LFEImpl >;
          for ( auto i : range( size() ) )
          {
            LFEImpl imp( i );
            map_[ i ].reset( new LFEObject( imp ) );
          }
        }
        /*
        else
        {
          DUNE_THROW(NotImplemented,"BDFMLocalFiniteElement not implemented for your choice." );
        }
        */
      }

      const LocalFiniteElementType& map( const std::size_t i ) const
      {
        assert( map_[ i ] );
        return *map_[ i ];
      }


      UniqueFacetOrientation< GridPartType > orientation_;
      std::array< std::unique_ptr< LocalFiniteElementType >, size() > map_;

      const int order_;
    };


    // BrezziDouglasFortinMariniSpace
    // ------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    using BrezziDouglasFortinMariniSpace
        = LocalFiniteElementSpace< BrezziDouglasFortinMariniLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >, FunctionSpace, Storage >;

    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    using DynamicBrezziDouglasFortinMariniSpace
        = LocalFiniteElementSpace< BrezziDouglasFortinMariniLocalFiniteElementMap< GridPart, FunctionSpace >, FunctionSpace, Storage >;


  } // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASFORTINMARINI_HH
