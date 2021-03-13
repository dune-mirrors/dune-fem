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
#include <dune/localfunctions/raviartthomas/raviartthomas0cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas0cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas12d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas1cube3d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas2cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas3cube2d.hh>
#include <dune/localfunctions/raviartthomas/raviartthomas4cube2d.hh>

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
        static_assert( AlwaysFalse< DomainField >::value, "RaviartThomasLocalFiniteElement not implemented for your choice." );
      };

      // 2d, Simplex, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::simplex( 2 ).id(), D, R, 2, 0 >
        : public RT02DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 8;
        using RT02DLocalFiniteElement< D, R >::RT02DLocalFiniteElement;
      };

      // 2d, Simplex, 1st order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::simplex( 2 ).id(), D, R, 2, 1 >
        : public RT12DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 8;
        using RT12DLocalFiniteElement< D, R >::RT12DLocalFiniteElement;
      };

      // 2d, Cube, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 0 >
        : public RT0Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16;
        using RT0Cube2DLocalFiniteElement< D, R >::RT0Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 1st order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 1 >
        : public RT1Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16;
        using RT1Cube2DLocalFiniteElement< D, R >::RT1Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 2nd order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 2 >
        : public RT2Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16;
        using RT2Cube2DLocalFiniteElement< D, R >::RT2Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 3rd order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 3 >
        : public RT3Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16;
        using RT3Cube2DLocalFiniteElement< D, R >::RT3Cube2DLocalFiniteElement;
      };

      // 2d, Cube, 4th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), D, R, 2, 4 >
        : public RT4Cube2DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 16;
        using RT4Cube2DLocalFiniteElement< D, R >::RT4Cube2DLocalFiniteElement;
      };

      // 3d, Cube, 0th order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 3 ).id(), D, R, 3, 0 >
        : public RT0Cube3DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 64;
        using RT0Cube3DLocalFiniteElement< D, R >::RT0Cube3DLocalFiniteElement;
      };

      // 3d, Cube, 1st order
      template< class D, class R >
      struct RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 3 ).id(), D, R, 3, 1 >
        : public RT1Cube3DLocalFiniteElement< D, R >
      {
        static constexpr std::size_t numOrientations = 64;
        using RT1Cube3DLocalFiniteElement< D, R >::RT1Cube3DLocalFiniteElement;
      };

    } // namespace Impl


    // RaviartThomasLocalFiniteElementMap
    // ----------------------------------

    template< class GridPart, class FunctionSpace, int polOrder >
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

      using KeyType = std::tuple<>;

      using TransformationType = PiolaTransformation< Geometry, FunctionSpace::dimRange >;

      using LocalFiniteElementType =
          Impl::RaviartThomasLocalFiniteElement< topologyId, DomainFieldType, RangeFieldType, dimLocal, polOrder >;

      using LocalBasisType          = typename LocalFiniteElementType::Traits::LocalBasisType;
      using LocalCoefficientsType   = typename LocalFiniteElementType::Traits::LocalCoefficientsType;
      using LocalInterpolationType  = typename LocalFiniteElementType::Traits::LocalInterpolationType;

      static constexpr auto size () -> std::size_t { return LocalFiniteElementType::numOrientations; }

      template< class ... Args >
      RaviartThomasLocalFiniteElementMap ( const GridPartType& gridPart, Args&& ... )
        : orientation_( gridPart )
      {
        for ( auto i : range( size() ) )
          map_[ i ] = LocalFiniteElementType( i );
      }

      int order () const { return polOrder; }

      template< class Entity >
      int order ( const Entity& entity ) const { return order; }

      template< class Entity >
      auto operator() ( const Entity& entity ) const
        -> std::tuple< std::size_t, const LocalBasisType&, const LocalInterpolationType& >
      {
        auto o = orientation_( entity );
        return std::tuple< std::size_t, const LocalBasisType&, const LocalInterpolationType& >(
            static_cast< std::size_t >( o ), map_[ o ].localBasis(), map_[ o ].localInterpolation() );
      }

      auto localCoefficients ( const GeometryType& type ) const
        -> const LocalCoefficientsType&
      {
        return map_[ 0 ].localCoefficients();
      }

      bool hasCoefficients ( const GeometryType& type ) const { return type == GeometryType( topologyId, dimLocal ); }
      auto gridPart () const -> const GridPartType& { return orientation_.gridPart(); }

    private:
      UniqueFacetOrientation< GridPartType > orientation_;
      std::array< LocalFiniteElementType, LocalFiniteElementType::numOrientations > map_;
    };


    // RaviartThomasSpace
    // ------------------

    template< class FunctionSpace, class GridPart, int polOrder, class Storage = CachingStorage >
    using RaviartThomasSpace
        = LocalFiniteElementSpace< RaviartThomasLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >, FunctionSpace, Storage >;


  } // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_RAVIARTTHOMAS_SPACE_HH
