#ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH
#define DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH

#include <cassert>
#include <cstddef>

#include <type_traits>
#include <utility>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // FiniteVolumeBasisFunctionSet
    // ----------------------------

    template< class Entity, class Range >
    struct FiniteVolumeBasisFunctionSet
    {
      /** \copydoc Dune::Fem::BasisFunctionSet::EntityType */
      typedef Entity EntityType;

      /** \copydoc Dune::Fem::BasisFunctionSet::FunctionSpaceType */
      typedef FunctionSpace< typename Entity::Geometry::ctype, typename Range::value_type,
                             Entity::Geometry::coorddimension, Range::dimension
                           > FunctionSpaceType;

      /** \copydoc Dune::Fem::BasisFunctionSet::DomainType */
      typedef typename FunctionSpaceType::DomainType DomainType;
      /** \copydoc Dune::Fem::BasisFunctionSet::RangeType */
      typedef typename FunctionSpaceType::RangeType RangeType;
      /** \copydoc Dune::Fem::BasisFunctionSet::JacobianRangeType */
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      /** \copydoc Dune::Fem::BasisFunctionSet::HessianRangeType */
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      /** \copydoc Dune::Fem::BasisFunctionSet::ReferenceElementType */
      typedef std::decay_t< decltype( Dune::ReferenceElements< typename EntityType::Geometry::ctype, EntityType::Geometry::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) ) > ReferenceElementType;

      /** \name Construction
       *  \{
       */

      FiniteVolumeBasisFunctionSet () : entity_( nullptr ) {}

      explicit FiniteVolumeBasisFunctionSet ( const EntityType &entity )
        : entity_( &entity )
      {}

      /** \} */

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::BasisFunctionSet::order */
      static constexpr int order () { return 0; }

      /** \copydoc Dune::Fem::BasisFunctionSet::size */
      static constexpr std::size_t size () { return RangeType::dimension; }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Quadrature, class Vector, class DofVector >
      void axpy ( const Quadrature &quadrature, const Vector &values, DofVector &dofs ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          axpy( quadrature[ qp ], values[ qp ], dofs );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
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

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          dofs[ i ] += valueFactor[ i ];
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {}

      /** \copydoc Dune::Fem::BasisFunctionSet::axpy */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::evaluateAll */
      template< class Quadrature, class DofVector, class RangeArray >
      void evaluateAll ( const Quadrature &quadrature, const DofVector &dofs, RangeArray &ranges ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          evaluateAll( quadrature[ qp ], dofs, ranges[ qp ] );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::evaluateAll */
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          value[ i ] = dofs[ i ];
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::evaluateAll */
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
        {
          values[ i ] = RangeType( 0 );
          values[ i ][ i ] = typename RangeType::field_type( 1 );
        }
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::jacobianAll */
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quadrature, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          jacobianAll( quadrature[ qp ], dofs, jacobians[ qp ] );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::jacobianAll */
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 0 );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::jacobianAll */
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          jacobians[ i ] = JacobianRangeType( 0 );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::hessianAll */
      template< class QuadratureType, class DofVector, class HessianArray >
      void hessianAll ( const QuadratureType &quadrature, const DofVector &dofs, HessianArray &hessians ) const
      {
        assert( hessians.size() >= quadrature.nop() );
        const unsigned int nop = quadrature.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          hessians[qp] = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::hessianAll */
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::hessianAll */
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        for( int i = 0; i < RangeType::dimension; ++i )
          hessians[ i ] = HessianRangeType( typename HessianRangeType::value_type( 0 ) );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::entity */
      const EntityType &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::referenceElementType */
      auto referenceElement () const
        -> decltype( Dune::ReferenceElements< typename EntityType::Geometry::ctype, EntityType::Geometry::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) )
      {
        return Dune::ReferenceElements< typename EntityType::Geometry::ctype, EntityType::Geometry::coorddimension >::general( type() );
      }

      /** \copydoc Dune::Fem::BasisFunctionSet::type */
      Dune::GeometryType type () const { return entity().type(); }

      //! \brief return true if entity pointer is set
      bool valid () const { return bool(entity_); }

    private:
      const EntityType *entity_;
    };

    /** \} */

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FINITEVOLUME_BASISFUNCTIONSET_HH
