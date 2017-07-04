#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALFEMAP_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_LOCALFEMAP_HH

// dune-localfunctions includes
#include <dune/localfunctions/rannacherturek.hh>

namespace Dune
{

  namespace Fem
  {

    template< class GridPart, class FunctionSpace >
    class RannacherTurekLocalFiniteElementMap
    {
      typedef RannacherTurekLocalFiniteElementMap< GridPart, FunctionSpace > ThisType;

    public:
      typedef GridPart GridPartType;
      typedef std::tuple< > KeyType;

      typedef typename FunctionSpace::DomainFieldType DomainFieldType;
      typedef typename FunctionSpace::RangeFieldType RangeFieldType;

      static const int dimLocal = GridPart::dimension;

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
      std::tuple< std::size_t, const LocalBasisType &, const LocalInterpolationType & > operator() ( const Entity &e ) const
      {
        return std::make_tuple(
          static_cast< std::size_t >( 0 ),
          localFe_.localBasis(),
          localFe_.localInterpolation() );
      }

      bool hasCoefficient ( const GeometryType &type ) const { return type.isCube(); }

      template< class GeometryType >
      LocalCoefficientsType localCoefficients ( const GeometryType &type ) const { return localFe_.localCoefficients(); }

      const GridPartType &gridPart () const { return gridPart_; }

    private:
      LocalFiniteElementType localFe_;
      const GridPartType &gridPart_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_LOCALFEMAP_HH
