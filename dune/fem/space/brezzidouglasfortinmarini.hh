#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASFORTINMARINI_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASFORTINMARINI_HH

#if HAVE_DUNE_LOCALFUNCTIONS

// C++ includes
#include <array>
#include <tuple>

// dune-common includes
#include <dune/common/power.hh>
#include <dune/common/typetraits.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-localfunction includes
#include <dune/localfunctions/brezzidouglasfortinmarini/bdfmcube.hh>

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

    template< class GridPart, class FunctionSpace, int polOrder >
    class BrezziDouglasFortinMariniLocalFiniteElementMap
    {
      using hasSingleGeometryType = GridPartCapabilities::hasSingleGeometryType< GridPart >;

      static_assert( hasSingleGeometryType::v, "`GridPart` has more the one geometry type." );

      static constexpr int dimLocal = GridPart::dimension;
      static constexpr unsigned int topologyId = hasSingleGeometryType::topologyId;
      static_assert( topologyId == Dune::Impl::CubeTopology< dimLocal >::type::id,
                     "Only defined for cube grids." );

      static_assert( dimLocal == FunctionSpace::dimRange, "`dimRange` has to be equal to `GridPart::dimension`" );

      using Geometry = typename GridPart::template Codim< 0 >::EntityType::Geometry;

    public:
      using GridPartType = GridPart;

      using DomainFieldType = typename FunctionSpace::DomainFieldType;
      using RangeFieldType = typename FunctionSpace::RangeFieldType;

      using KeyType = std::tuple<>;

      using TransformationType = PiolaTransformation< Geometry, FunctionSpace::dimRange >;

      using LocalFiniteElementType = BDFMCubeLocalFiniteElement< DomainFieldType, RangeFieldType, dimLocal, polOrder >;

      using LocalBasisType          = typename LocalFiniteElementType::Traits::LocalBasisType;
      using LocalCoefficientsType   = typename LocalFiniteElementType::Traits::LocalCoefficientsType;
      using LocalInterpolationType  = typename LocalFiniteElementType::Traits::LocalInterpolationType;

      static constexpr auto size () -> std::size_t { return StaticPower< 2, 2*dimLocal >::power; }

      template< class ... Args >
      BrezziDouglasFortinMariniLocalFiniteElementMap ( const GridPartType& gridPart, Args&& ... )
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
      std::array< LocalFiniteElementType, size() > map_;
    };


    // BrezziDouglasFortinMariniSpace
    // ------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage = CachingStorage >
    using BrezziDouglasFortinMariniSpace
        = LocalFiniteElementSpace< BrezziDouglasFortinMariniLocalFiniteElementMap< GridPart, FunctionSpace, polOrder >, FunctionSpace, Storage >;


  } // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_LOCALFUNCTIONS

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASFORTINMARINI_HH
