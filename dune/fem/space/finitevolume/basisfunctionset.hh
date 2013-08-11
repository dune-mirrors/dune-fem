#ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH

#include <cassert>
#include <cstddef>

#include <dune/common/nullptr.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeBasisFunctionSet
    // ----------------------------

    template< class Entity, class Range >
    struct FiniteVolumeBasisFunctionSet
    {
      typedef Entity EntityType;

      typedef FunctionSpace< typename Entity::Geometry::ctype, typename Range::value_type, 
                             Entity::Geometry::coorddimension, Range::dimension
                           > FunctionSpaceType;
       
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef Dune::ReferenceElement< typename DomainType::value_type, DomainType::dimension > ReferenceElementType;

      FiniteVolumeBasisFunctionSet () : entity_( nullptr ) {}

      FiniteVolumeBasisFunctionSet ( const EntityType &entity ) : entity_( &entity ) {}

      int order () const { return 0; }

      std::size_t size () const { return RangeType::dimension; }

      template< class Quadrature, class Vector, class DofVector >
      void axpy ( const Quadrature &quadrature, const Vector &values, DofVector &dofs ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          axpy( quadrature[ qp ], values[ qp ], dofs );
      }

      template< class Quadrature, class VectorA, class VectorB, class DofVector >
      void axpy ( const Quadrature &quadrature, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          axpy( quadrature[ qp ], valuesA[ qp ], dofs );
          axpy( quadrature[ qp ], valuesB[ qp ], dofs );
        }
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          dofs[ i ] += valueFactor[ i ];
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {}

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
      }

      template< class Quadrature, class DofVector, class RangeArray >
      void evaluateAll ( const Quadrature &quadrature, const DofVector &dofs, RangeArray &ranges ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          evaluateAll( quadrature[ qp ], dofs, ranges[ qp ] );
      }

      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          value[ i ] = dofs[ i ];
      }

      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
        {
          values[ i ] = RangeType( 0 );
          values[ i ][ i ] = typename RangeType::field_type( 1 );
        }
      }

      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quadrature, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          jacobianAll( quadrature[ qp ], dofs, jacobians[ qp ] );
      }

      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 0 );
      }

      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          jacobians[ i ] = JacobianRangeType( 0 );
      }

      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          hessians[ i ] = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      const ReferenceElementType &referenceElement () const
      {
        return Dune::ReferenceElements< typename EntityType::ctype, EntityType::dimensionworld >::general( type() );
      }

      Dune::GeometryType type () const { return entity().type(); }

    private:
      const EntityType *entity_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH
