#ifndef DUNE_FEM_SPACE_RAVIARTTHOMAS_SPACE_HH
#define DUNE_FEM_SPACE_RAVIARTTHOMAS_SPACE_HH

#if HAVE_DUNE_LOCALFUNCTIONS

// C++ includes
#include <array>
#include <tuple>

// dune-common includes
#include <dune/common/typetraits.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-localfunction includes
#include <dune/localfunctions/raviartthomas/raviartthomas02d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas03d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas2cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas3cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas4cube2d.hh>
#include <dune/localfunctions/common/virtualwrappers.hh>

// dune-fem includes
#include <dune/fem/space/common/uniquefacetorientation.hh>
#include <dune/fem/space/basisfunctionset/piolatransformation.hh>
#include <dune/fem/space/localfiniteelement/space.hh>


namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {
      // RaviartThomasLocalFiniteElement
      // -------------------------------

      template< unsigned int id, class DomainField, class RangeField, int dimension, int order >
      struct RaviartThomasLocalFiniteElement
      {
        // 0 here means that the element is not implemented
        static constexpr std::size_t numOrientations = 0;

        RaviartThomasLocalFiniteElement()
        {
          DUNE_THROW(NotImplemented,"RaviartThomasLocalFiniteElement not implemented for your choice." );
        }
      };

      // 2d, Simplex, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::simplex( 2 ).id(), D, R, 2, 0 >
        : public RT02DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 8; // 2^3
        using RT02DLocalFiniteElement< D, R >::RT02DLocalFiniteElement;
      };

      // 2d, Simplex, 1st order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::simplex( 2 ).id(), D, R, 2, 1 >
        : public RT12DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 8; // 2^3
        using RT12DLocalFiniteElement< D, R >::RT12DLocalFiniteElement;
      };

      // 3d, Simplex, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::simplex( 3 ).id(), D, R, 3, 0 >
        : public RT03DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16; // 2^4
        using RT03DLocalFiniteElement< D, R >::RT03DLocalFiniteElement;
      };

      // 2d, Cube, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 0 >
        : public RT0Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16; // 2^4
        using RT0Cube2DLocalFiniteElement< D, R >::RT0Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 1st order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 1 >
        : public RT1Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16; // 2^4
        using RT1Cube2DLocalFiniteElement< D, R >::RT1Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 2nd order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 2 >
        : public RT2Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16; // 2^4
        using RT2Cube2DLocalFiniteElement< D, R >::RT2Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 3rd order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 3 >
        : public RT3Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16; // 2^4
        using RT3Cube2DLocalFiniteElement< D, R >::RT3Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 4th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 4 >
        : public RT4Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16; // 2^4
        using RT4Cube2DLocalFiniteElement< D, R >::RT4Cube2DLocalFiniteElement;
      };

      // 3d, Cube, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 3 ).id(), D, R, 3, 0 >
        : public RT0Cube3DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 64; // 2^6
        using RT0Cube3DLocalFiniteElement< D, R >::RT0Cube3DLocalFiniteElement;
      };

      // 3d, Cube, 1st order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 3 ).id(), D, R, 3, 1 >
        : public RT1Cube3DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 64; // 2^6
        using RT1Cube3DLocalFiniteElement< D, R >::RT1Cube3DLocalFiniteElement;
      };

    } // namespace Impl

    // RaviartThomasLocalFiniteElementMap
    // ----------------------------------

    template< class GridPart, class FunctionSpace, int polOrder = -1 >
    class RaviartThomasLocalFiniteElementMap
    {
      using hasSingleGeometryType = GridPartCapabilities::hasSingleGeometryType< GridPart >;

      static_assert( hasSingleGeometryType::v, "`GridPart` has more the one geometry type." );

      static constexpr int dimLocal = GridPart::dimension;
      static constexpr unsigned int topologyId = hasSingleGeometryType::topologyId;

      static_assert( dimLocal == FunctionSpace::dimRange, "`dimRange` has to be equal to `GridPart::dimension`" );

      using Geometry = typename GridPart::template Codim< 0 >::EntityType::Geometry;

    public:
      using GridPartType = GridPart;

      using DomainFieldType = typename FunctionSpace::DomainFieldType;
      using RangeFieldType = typename FunctionSpace::RangeFieldType;

      typedef unsigned int KeyType;

      using TransformationType = PiolaTransformation< Geometry, FunctionSpace::dimRange >;

    protected:
      // for dynamic space the polOrder is passed as negative value
      static constexpr bool dynamicOrder = polOrder < 0;

      using DomainType = typename FunctionSpace::DomainType;
      using RangeType  = typename FunctionSpace::RangeType;
      typedef LocalBasisTraits<DomainFieldType, FunctionSpace::dimDomain, DomainType,
                               RangeFieldType,  FunctionSpace::dimRange,  RangeType,
                               typename FunctionSpace::JacobianRangeType > LBT;

      // this is the same for each order, so take order 0 which exists for all combinations
      static constexpr std::size_t numOrientations =
          Impl::RaviartThomasLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, 0 >::numOrientations;

    public:
      // type of LocalFiniteElement
      typedef typename std::conditional< dynamicOrder,
              // virtual LFE type
              LocalFiniteElementVirtualInterface< LBT >,
              // real implementation type
              Impl::RaviartThomasLocalFiniteElement< topologyId, DomainFieldType,
                  RangeFieldType, dimLocal, polOrder > > :: type     LocalFiniteElementType;

      using LocalBasisType          = typename LocalFiniteElementType::Traits::LocalBasisType;
      using LocalCoefficientsType   = typename LocalFiniteElementType::Traits::LocalCoefficientsType;
      using LocalInterpolationType  = typename LocalFiniteElementType::Traits::LocalInterpolationType;

      static constexpr auto size () -> std::size_t { return numOrientations; }

      RaviartThomasLocalFiniteElementMap ( const GridPartType& gridPart, const unsigned int ord )
        : orientation_( gridPart ), order_( dynamicOrder ? ord : polOrder )
      {
        if constexpr ( dynamicOrder )
        {
          if ( order_ == 0 )
            this->template createLFE< 0 > ();
          else if ( order_ == 1 )
            this->template createLFE< 1 > ();
          else if ( order_ == 2 )
            this->template createLFE< 2 > ();
          else if ( order_ == 3 )
            this->template createLFE< 3 > ();
          else if ( order_ == 4 )
            this->template createLFE< 4 > ();
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

      auto localCoefficients ( const GeometryType& type ) const -> const LocalCoefficientsType&
      {
        return map( 0 ).localCoefficients();
      }

      bool hasCoefficients ( const GeometryType& type ) const { return type == GeometryType( topologyId, dimLocal ); }
      auto gridPart () const -> const GridPartType& { return orientation_.gridPart(); }

    private:
      template <int p>
      void createLFE()
      {
        using LFEImpl = Impl::RaviartThomasLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, p >;
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
          DUNE_THROW(NotImplemented,"RaviartThomasLocalFiniteElement not implemented for your choice." );
        }
      }

      const LocalFiniteElementType& map( const size_t i ) const
      {
        assert( map_[i] );
        return *(map_[i]);
      }

      UniqueFacetOrientation< GridPartType > orientation_;
      std::array< std::unique_ptr< LocalFiniteElementType >, numOrientations > map_;

      const int order_;
    };


    // RaviartThomasSpace
    // ------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    using RaviartThomasSpace
        = LocalFiniteElementSpace< RaviartThomasLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >, FunctionSpace, Storage >;

    // dynamic polynomial order choice in this case
    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    using DynamicRaviartThomasSpace
        = LocalFiniteElementSpace< RaviartThomasLocalFiniteElementMap< GridPart, FunctionSpace >, FunctionSpace, Storage >;


  } // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_RAVIARTTHOMAS_SPACE_HH
