#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_HH

#if HAVE_DUNE_LOCALFUNCTIONS
// dune-localfunctions includes

#include <dune/localfunctions/rannacherturek.hh>
#include <dune/fem/gridpart/common/capabilities.hh>

#include <dune/fem/space/localfiniteelement/space.hh>

#include <dune/fem/space/rannacherturek/space.hh>

namespace Dune
{
  namespace Fem
  {

    template< class FunctionSpace, class GridPart >
    class RannacherTurekLocalFiniteElementMap
    {
      static const int dimLocal = GridPart::dimension;

      typedef RannacherTurekLocalFiniteElementMap< FunctionSpace, GridPart > ThisType;
      static_assert( Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::v,
                     "GridPart has more than one geometry type." );
      static_assert( ((Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPart >::topologyId
                     ^ ((1 << dimLocal)-1)) >> 1 == 0),
                     "Space only defined on cube grids." );
    public:
      typedef GridPart GridPartType;

      typedef std::tuple< > KeyType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      typedef RannacherTurekLocalFiniteElement< DomainFieldType, RangeFieldType, dimLocal > LocalFiniteElementType;
      typedef typename LocalFiniteElementType::Traits::LocalBasisType LocalBasisType;
      typedef typename LocalFiniteElementType::Traits::LocalCoefficientsType LocalCoefficientsType;
      typedef typename LocalFiniteElementType::Traits::LocalInterpolationType LocalInterpolationType;

      template< class ... Args >
      RannacherTurekLocalFiniteElementMap ( const GridPart &gridPart, Args ... args )
        : gridPart_( gridPart ) {}

      static std::size_t size () { return 1; }

      int order () const { return localFe_.localBasis().order(); }

      template< class Entity >
      int order ( const Entity &entity ) const { return order(); }

      template< class Entity >
      auto operator() ( const Entity &e ) const
      {
        return std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & >
        { static_cast< std::size_t >( 0 ),
          localFe_.localBasis(),
          localFe_.localInterpolation() };
      }

      bool hasCoefficients ( const GeometryType &type ) const { return type.isCube(); }

      const LocalCoefficientsType& localCoefficients ( const GeometryType &type ) const
      {
        return localFe_.localCoefficients();
      }

      const GridPartType &gridPart () const { return gridPart_; }

    private:
      LocalFiniteElementType localFe_;
      const GridPartType &gridPart_;
    };

    template< class FunctionSpace, class GridPart, class Storage = CachingStorage >
    using RannacherTurekSpace
    = LocalFiniteElementSpace<
    RannacherTurekLocalFiniteElementMap< FunctionSpace, GridPart >, FunctionSpace, Storage >;

  }
}

#endif // HAVE_DUNE_LOCALFUNCTIONS
#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_HH
