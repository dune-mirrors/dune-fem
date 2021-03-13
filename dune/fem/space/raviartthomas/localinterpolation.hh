#ifndef DUNE_FEM_SPACE_RAVIARTTHOMAS_LOCALINTERPOLATION_HH
#define DUNE_FEM_SPACE_RAVIARTTHOMAS_LOCALINTERPOLATION_HH

// alternative interpolation used for testing

// C++ includes
#include <cassert>
#include <vector>
#include <utility>

// dune-common includes
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/quadraturerules.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/piolatransformation.hh>


namespace Dune
{

  namespace Fem
  {

    namespace Impl
    {

      // forward declarations
      // --------------------

      template< unsigned int, class, class, int, int > struct RaviartThomasLocalFiniteElement;


      // RaviartThomasLocalInterpolationBasis
      // ------------------------------------

      /*
       *  These are mostly copies from the interpolation implementations in dune-localfunctions
       *  (dune/localfunctions/raviartthomas/raviartthomas * /raviartthomas * localinterpolation.hh)
       */

      template< class LocalFiniteElement >
      struct RaviartThomasLocalInterpolationBasis
      {
        static_assert( AlwaysFalse< LocalFiniteElement >::value,
                       "`RaviartThomasLocalInterpolationBasis` not implemented for your choice." );
      };


      // 0th order
      template< unsigned int id, class Domain, class Range, int dim >
      struct RaviartThomasLocalInterpolationBasis< RaviartThomasLocalFiniteElement< id, Domain, Range, dim, 0 > >
      {
        using DomainType = FieldVector< Domain, dim >;
        using FaceDomainType = FieldVector< Domain, dim-1 >;
        using RangeType = FieldVector< Range, dim >;
        using RangeFieldType = Range;

        RaviartThomasLocalInterpolationBasis () = default;
        explicit RaviartThomasLocalInterpolationBasis ( unsigned int orientations ) : orientations_( orientations ) {}

        void trace ( int facet, const FaceDomainType& x, std::vector< std::pair< int, RangeFieldType > >& basis ) const
        {
          assert( basis.size() >= size( 1 ) );
          basis[ 0 ] = std::make_pair( facet, sign( facet ) );
        }

        void interior ( const DomainType& x, std::vector< std::pair< int, RangeType > >& basis ) const {}

        constexpr auto size ( int codim ) const -> std::size_t { assert( codim < 2 ); return (codim == 0) ? 0 : 1; }
        constexpr int order ( int codim ) const { assert( codim < 2 ); return (codim == 0) ? 0 : 1; }

      private:
        auto sign ( int facet ) const -> RangeFieldType { return (orientations_ & (1u << facet)) ? -1.0 : 1.0; }

        unsigned int orientations_ = 0;
      };

      // 1st order, 2d Simplex
      template< class Domain, class Range >
      struct RaviartThomasLocalInterpolationBasis< RaviartThomasLocalFiniteElement< Dune::Impl::SimplexTopology< 2 >::type::id, Domain, Range, 2, 1 > >
      {
        using DomainType = FieldVector< Domain, 2 >;
        using FaceDomainType = FieldVector< Domain, 1 >;
        using RangeType = FieldVector< Range, 2 >;
        using RangeFieldType = Range;

        RaviartThomasLocalInterpolationBasis () = default;
        explicit RaviartThomasLocalInterpolationBasis ( unsigned int orientations ) : orientations_( orientations ) {}

        void trace ( int facet, const FaceDomainType& x, std::vector< std::pair< int, RangeFieldType > >& basis ) const
        {
          assert( basis.size() >= size( 1 ) );
          const RangeFieldType temp = (facet==1) ? 1.0 - 2.0*x : 2.0*x - 1.0;
          basis[ 0 ] = std::make_pair( facet, sign( facet ) );
          basis[ 1 ] = std::make_pair( facet+3, temp );
        }

        void interior ( const DomainType& x, std::vector< std::pair< int, RangeType > >& basis ) const
        {
          assert( basis.size() >= size( 0 ) );
          basis[ 0 ] = std::make_pair( 6, RangeType{ 1.0, 0.0 } );
          basis[ 1 ] = std::make_pair( 7, RangeType{ 0.0, 1.0 } );
        }

        constexpr auto size ( int codim ) const -> std::size_t { assert( codim < 2 ); return (codim == 0) ? 2 : 2; }
        constexpr int order ( int codim ) const { assert( codim < 2 ); return (codim == 0) ? 8 : 4; }


      private:
        auto sign ( int facet ) const -> RangeFieldType { return (orientations_ & (1u << facet)) ? -1.0 : 1.0; }

        unsigned int orientations_ = 0;
      };

      // 1st order, 2d Cube
      template< class Domain, class Range >
      struct RaviartThomasLocalInterpolationBasis< RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 2 ).id(), Domain, Range, 2, 1 > >
      {
        using DomainType = FieldVector< Domain, 2 >;
        using FaceDomainType = FieldVector< Domain, 1 >;
        using RangeType = FieldVector< Range, 2 >;
        using RangeFieldType = Range;

        RaviartThomasLocalInterpolationBasis () = default;
        explicit RaviartThomasLocalInterpolationBasis ( unsigned int orientations ) : orientations_( orientations ) {}

        void trace ( int facet, const FaceDomainType& x, std::vector< std::pair< int, RangeFieldType > >& basis ) const
        {
          assert( basis.size() >= size( 1 ) );
          const RangeFieldType temp = (facet > 0 && facet < 3 ) ? 1.0 - 2.0*x : 2.0*x - 1.0;
          basis[ 0 ] = std::make_pair( 2*facet  , sign( facet ) );
          basis[ 1 ] = std::make_pair( 2*facet+1, temp );
        }

        void interior ( const DomainType& x, std::vector< std::pair< int, RangeType > >& basis ) const
        {
          assert( basis.size() >= size( 0 ) );
          basis[ 0 ] = std::make_pair(  8, RangeType{  1.0,  0.0 } );
          basis[ 1 ] = std::make_pair(  9, RangeType{  0.0,  1.0 } );
          basis[ 2 ] = std::make_pair( 10, RangeType{ x[1],  0.0 } );
          basis[ 3 ] = std::make_pair( 11, RangeType{  0.0, x[0] } );
        }

        constexpr auto size ( int codim ) const -> std::size_t { assert( codim < 2 ); return (codim == 0) ? 4 : 2; }
        constexpr int order ( int codim ) const { assert( codim < 2 ); return (codim == 0) ? 3 : 3; }

      private:
        auto sign ( int facet ) const -> RangeFieldType { return (orientations_ & (1u << facet)) ? -1.0 : 1.0; }

        unsigned int orientations_ = 0;
      };

      // 1st order, 3d Cube
      template< class Domain, class Range >
      struct RaviartThomasLocalInterpolationBasis< RaviartThomasLocalFiniteElement< Dune::GeometryTypes::cube( 3 ).id(), Domain, Range, 3, 1 > >
      {
        using DomainType = FieldVector< Domain, 3 >;
        using FaceDomainType = FieldVector< Domain, 2 >;
        using RangeType = FieldVector< Range, 3 >;
        using RangeFieldType = Range;

        RaviartThomasLocalInterpolationBasis () = default;
        explicit RaviartThomasLocalInterpolationBasis ( unsigned int orientations ) : orientations_( orientations ) {}

        void trace ( int facet, const FaceDomainType& x, std::vector< std::pair< int, RangeFieldType > >& basis ) const
        {
          assert( basis.size() >= size( 1 ) );

          const RangeFieldType tempX = 2.0*x[0] - 1.0;
          const RangeFieldType tempY = 2.0*x[1] - 1.0;

          basis[ 0 ] = std::make_pair( facet, sign( facet ) );

          switch( facet )
          {
            case 0:
            case 5:
              basis[ 1 ] = std::make_pair( facet+ 6,        tempX );
              basis[ 2 ] = std::make_pair( facet+12,        tempY );
              basis[ 3 ] = std::make_pair( facet+18,  tempX*tempY );
              break;
            case 1:
            case 4:
              basis[ 1 ] = std::make_pair( facet+ 6,       -tempX );
              basis[ 2 ] = std::make_pair( facet+12,       -tempY );
              basis[ 3 ] = std::make_pair( facet+18, -tempX*tempY );
              break;
            case 2:
              basis[ 1 ] = std::make_pair( facet+ 6,       -tempX );
              basis[ 2 ] = std::make_pair( facet+12,        tempY );
              basis[ 3 ] = std::make_pair( facet+18, -tempX*tempY );
              break;
            case 3:
              basis[ 1 ] = std::make_pair( facet+ 6,        tempX );
              basis[ 2 ] = std::make_pair( facet+12,       -tempY );
              basis[ 3 ] = std::make_pair( facet+18,  tempX*tempY );
              break;
          }
        }

        void interior ( const DomainType& x, std::vector< std::pair< int, RangeType > >& basis ) const
        {
          assert( basis.size() >= size( 0 ) );
          basis[  0 ] = std::make_pair( 24, RangeType{       1.0,       0.0,       0.0 } );
          basis[  1 ] = std::make_pair( 25, RangeType{       0.0,       1.0,       0.0 } );
          basis[  2 ] = std::make_pair( 26, RangeType{       0.0,       0.0,       1.0 } );
          basis[  3 ] = std::make_pair( 27, RangeType{      x[1],       0.0,       0.0 } );
          basis[  4 ] = std::make_pair( 28, RangeType{      x[2],       0.0,       0.0 } );
          basis[  5 ] = std::make_pair( 29, RangeType{       0.0,      x[0],       0.0 } );
          basis[  6 ] = std::make_pair( 30, RangeType{       0.0,      x[2],       0.0 } );
          basis[  7 ] = std::make_pair( 31, RangeType{       0.0,       0.0,      x[0] } );
          basis[  8 ] = std::make_pair( 32, RangeType{       0.0,       0.0,      x[1] } );
          basis[  9 ] = std::make_pair( 33, RangeType{ x[1]*x[2],       0.0,       0.0 } );
          basis[ 10 ] = std::make_pair( 34, RangeType{       0.0, x[0]*x[2],       0.0 } );
          basis[ 11 ] = std::make_pair( 35, RangeType{       0.0,       0.0, x[0]*x[1] } );
        }

        constexpr auto size ( int codim ) const -> std::size_t { assert( codim < 2 ); return (codim == 0) ? 12 : 4; }
        constexpr int order ( int codim ) const { assert( codim < 2 ); return (codim == 0) ? 3 : 3; }

      private:
        auto sign ( int facet ) const -> RangeFieldType { return (orientations_ & (1u << facet)) ? -1.0 : 1.0; }

        unsigned int orientations_ = 0;
      };

    }

    template< class GridPart, class LocalFiniteElement, int dimRange = GridPart::dimension >
    class RaviartThomasLocalInterpolation
    {
    public:
      using GridPartType = GridPart;
      using LocalFiniteElementType = LocalFiniteElement;

      using EntityType = typename GridPartType::template Codim< 0 >::EntityType;

    protected:
      using LocalInterpolationBasisType = Impl::RaviartThomasLocalInterpolationBasis< LocalFiniteElementType >;

      using Geometry = typename EntityType::Geometry;
      using ReferenceElementType = ReferenceElement< Geometry >;

      using TransformationType = PiolaTransformation< Geometry, dimRange >;
      using InverseTransformationType = typename TransformationType::InverseTransformationType;

      using RangeType = typename LocalInterpolationBasisType::RangeType;
      using RangeFieldType = typename LocalInterpolationBasisType::RangeFieldType;

      using VolumeQuadratures = QuadratureRules< typename Geometry::ctype, ReferenceElementType::dimension >;
      using FaceQuadratures = QuadratureRules< RangeFieldType, ReferenceElementType::dimension-1 >;

    public:
      explicit RaviartThomasLocalInterpolation ( const EntityType& entity )
        : RaviartThomasLocalInterpolation( entity, 0, -1 )
      {}

      RaviartThomasLocalInterpolation ( const EntityType& entity, unsigned int orientations, int order = -1 )
        : geometry_( entity.geometry() ),
          refElement_( Dune::referenceElement( geometry_ ) ),
          localBasis_( orientations ),
          order_( order )
      {}

      template< class LocalFunction, class LocalDofVector >
      void interior ( const LocalFunction& lf, LocalDofVector& dofs ) const
      {
        if ( !hasInterior() )
          return;

        assert( dofs.size() >= localBasis().size( 0 ) + referenceElement().size( 1 ) * localBasis().size( 1 ) );
        interior_.resize( localBasis().size( 0 ) );

        for( const auto& qp : getQuadrature( (order_ == -1) ? localBasis().order( 0 ) : order_ ) )
        {
          auto invPiola = inverseTransformation( qp.position() );

          localBasis().interior( qp.position(), interior_ );
          auto value = invPiola.apply( lf( qp.position() ) );

          for ( const auto& base : interior_ )
            dofs[ base.first ] += qp.weight() * (value * base.second);
        }
      }

      template< class LocalFunction, class LocalDofVector >
      void trace ( int facet, const LocalFunction& lf, LocalDofVector& dofs ) const
      {
        assert( dofs.size() >= localBasis().size( 0 ) + referenceElement().size( 1 ) * localBasis().size( 1 ) );
        traces_.resize( localBasis().size( 1 ) );

        auto embedding = referenceElement().template geometry< 1 >( facet );
        auto normal = referenceElement().integrationOuterNormal( facet );

        for ( const auto& qp : getQuadrature( facet, (order_ == -1) ? localBasis().order( 1 ) : order_ ) )
        {
          auto p = embedding.global( qp.position() );
          auto invPiola = inverseTransformation( p );

          localBasis().trace( facet, qp.position(), traces_ );
          auto value = invPiola.apply( lf( p ) );

          for ( const auto& base : traces_ )
            dofs[ base.first ] += qp.weight() * (value * normal) * base.second;
        }
      }

      template< class LocalFunction, class LocalDofVector >
      void interiorTrace ( int facet, const LocalFunction& lf, LocalDofVector& dofs ) const
      {
        if ( !hasInterior() )
          return;

        assert( dofs.size() >= localBasis().size( 0 ) + referenceElement().size( 1 ) * localBasis().size( 1 ) );
        interior_.resize( localBasis().size( 0 ) );

        auto embedding = referenceElement().template geometry< 1 >( facet );

        for( const auto& qp : getQuadrature( facet, (order_ == -1) ? localBasis().order( 0 ) : order_ ) )
        {
          auto p = embedding.global( qp.position() );
          auto invPiola = inverseTransformation( p );

          localBasis().interior( p, interior_ );
          auto value = invPiola.apply( lf( p ) );

          for ( const auto& base : interior_ )
            dofs[ base.first ] += qp.weight() * (value * base.second);
        }
      }

      template< class LocalFunction, class LocalDofVector >
      void operator() ( const LocalFunction& lf, LocalDofVector& dofs ) const
      {
        for ( int facet : range( referenceElement().size( 1 ) ) )
          trace( facet, lf, dofs );
        interior( lf, dofs );
      }

      bool hasInterior () const { return localBasis().size( 0 ) > 0; }

    protected:
      auto geometry () const -> const Geometry& { return geometry_; }
      auto referenceElement () const -> const ReferenceElementType& { return refElement_; }
      auto localBasis () const -> const LocalInterpolationBasisType& { return localBasis_; }

      template< class Point >
      auto transformation ( const Point& p ) const { return TransformationType( geometry(), p ); }

      template< class Point >
      auto inverseTransformation ( const Point& p ) const { return InverseTransformationType( geometry(), p ); }

      auto getQuadrature ( int order ) const { return VolumeQuadratures::rule( referenceElement().type(), order ); }
      auto getQuadrature ( int facet, int order ) const { return FaceQuadratures::rule( referenceElement().type( facet, 1 ), order ); }

    private:
      Geometry geometry_;
      ReferenceElementType refElement_;
      LocalInterpolationBasisType localBasis_;
      const int order_;

      mutable std::vector< std::pair< int, RangeFieldType > > traces_;
      mutable std::vector< std::pair< int, RangeType  > > interior_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RAVIARTTHOMAS_LOCALINTERPOLATION_HH
