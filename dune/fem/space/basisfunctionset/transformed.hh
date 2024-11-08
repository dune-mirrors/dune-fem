#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TRANSFORMED_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_TRANSFORMED_HH

// C++ includes
#include <cassert>
#include <cstddef>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/space/basisfunctionset/default.hh>
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/basisfunctionset/transformation.hh>
#include <dune/fem/space/common/functionspace.hh>


namespace Dune
{

  namespace Fem
  {

    // TransformedBasisFunctionSet
    // ---------------------------

    /**
     * \brief implementation of a basis function set for given entity
     *
     * \tparam  Entity            entity type
     * \tparam  ShapeFunctionSet  shape function set
     * \tparam  Transformation    transformation
     *
     * \note ShapeFunctionSet must be a copyable object. For most
     *       non-trivial implementations, you may want to use a
     *       proxy, see file
       \code
       <dune/fem/space/shapefunctionset/proxy.hh>
       \endcode
     */
    template< class Entity, class ShapeFunctionSet, class Transformation >
    class TransformedBasisFunctionSet
      : public BasisFunctionSetStorage< Entity, ShapeFunctionSet >
    {
      typedef BasisFunctionSetStorage< Entity, ShapeFunctionSet > BaseType;
      typedef TransformedBasisFunctionSet< Entity, ShapeFunctionSet, Transformation > ThisType;

    public:
      //! \brief entity type
      typedef typename BaseType :: EntityType  EntityType;

      //! \brief geometry
      typedef typename BaseType :: Geometry    Geometry;

      //! \brief shape function set type
      typedef typename BaseType :: ShapeFunctionSetType  ShapeFunctionSetType;

    protected:
      typedef Geometry GeometryType;
      typedef typename GeometryType::ctype ctype;
      typedef typename GeometryType::JacobianTransposed JacobianTransposed;

    public:
      //! \brief type of function space
      typedef typename ShapeFunctionSetType::FunctionSpaceType FunctionSpaceType;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      //! \brief type of reference element
      typedef typename BaseType :: ReferenceElementType  ReferenceElementType;

      //! \brief constructor
      TransformedBasisFunctionSet () {}

      //! \brief constructor
      explicit TransformedBasisFunctionSet ( const EntityType &entity, const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
        : BaseType( entity, shapeFunctionSet )
      {
      }

      TransformedBasisFunctionSet ( const TransformedBasisFunctionSet &other ) = default;
      TransformedBasisFunctionSet &operator= ( const TransformedBasisFunctionSet &other ) = default;

      using BaseType :: entity;
      using BaseType :: geometry;
      using BaseType :: type;
      using BaseType :: shapeFunctionSet;
      using BaseType :: order;
      using BaseType :: size;
      using BaseType :: referenceElement;

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class QuadratureType, class Vector, class DofVector >
      void axpy ( const QuadratureType &quad, const Vector &values, DofVector &dofs ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          axpy( quad[ qp ], values[ qp ], dofs );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       *
       *  \note valuesA and valuesB can be vectors of RangeType or
       *        JacobianRangeType
       */
      template< class QuadratureType, class VectorA, class VectorB, class DofVector >
      void axpy ( const QuadratureType &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          axpy( quad[ qp ], valuesA[ qp ], dofs );
          axpy( quad[ qp ], valuesB[ qp ], dofs );
        }
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        RangeType factor = transformation( coordinate( x ) ).apply_t( valueFactor );
        FunctionalAxpyFunctor< RangeType, DofVector > f( factor, dofs );
        shapeFunctionSet().evaluateEach( x, f );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        typedef typename GeometryType::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const Geometry &geo = geometry();
        const GeometryJacobianInverseTransposedType &gjit = geo.jacobianInverseTransposed( coordinate( x ) );
        JacobianRangeType tmpJacobianFactor( RangeFieldType( 0 ) );
        for( int r = 0; r < FunctionSpaceType::dimRange; ++r )
          gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

        tmpJacobianFactor = transformation( coordinate( x ) ).apply_t( tmpJacobianFactor );
        FunctionalAxpyFunctor< JacobianRangeType, DofVector > f( tmpJacobianFactor, dofs );
        shapeFunctionSet().jacobianEach( x, f );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy( const Point &x, const HessianRangeType &hessianFactor, DofVector &dofs ) const
      {
        DUNE_THROW( NotImplemented, "hessian axpy for TransformedBasisFunctionSet not implemented." );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      /** \copydoc BasisFunctionSet::evaluateAll( quad, dofs, ranges ) */
      template< class QuadratureType, class DofVector, class RangeArray >
      void evaluateAll ( const QuadratureType &quad, const DofVector &dofs, RangeArray &ranges ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        value = 0;
        AxpyFunctor< DofVector, RangeType > f( dofs, value );
        shapeFunctionSet().evaluateEach( x, f );
        value = transformation( coordinate( x ) ).apply( value );
      }

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        assert( values.size() >= size() );
        AssignFunctor< RangeArray > f( values );
        shapeFunctionSet().evaluateEach( x, f );
        auto transform = transformation( coordinate( x ) );
        for( auto &e : values )
          e = transform.apply( e );
      }

      /** \copydoc BasisFunctionSet::jacobianAll( quad, dofs, jacobians ) */
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        JacobianRangeType localJacobian( RangeFieldType( 0 ) );
        AxpyFunctor< DofVector, JacobianRangeType > f( dofs, localJacobian );
        shapeFunctionSet().jacobianEach( x, f );
        const Geometry &geo = geometry();
        JacobianTransformation< GeometryType >( geo, coordinate( x ) )( localJacobian, jacobian );
        jacobian = transformation( coordinate( x ) ).apply( jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        assert( jacobians.size() >= size() );
        const Geometry &geo = geometry();
        typedef JacobianTransformation< GeometryType > JacobianTransformation;
        JacobianTransformation jacobianTransformation( geo, coordinate( x ) );
        AssignFunctor< JacobianRangeArray, JacobianTransformation > f( jacobians, jacobianTransformation );
        shapeFunctionSet().jacobianEach( x, f );
        auto transform = transformation( coordinate( x ) );
        for( auto &jacobian : jacobians )
          jacobian = transform.apply( jacobian );
      }

      //! \todo please doc me
      template< class Point, class DofVector, class HessianRange >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRange& hessian ) const
      {
        if( ! std::is_same< HessianRange, HessianRangeType >:: value )
          DUNE_THROW( NotImplemented, "hessianAll: HessianRange mismatch!" );

        DUNE_THROW( NotImplemented, "hessianAll for TransformedBasisFunctionSet not implemented." );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        DUNE_THROW( NotImplemented, "hessianAll for TransformedBasisFunctionSet not implemented." );
      }

      // Non-interface methods
      // ---------------------

      Transformation transformation ( const DomainType& x ) const
      {
        return Transformation( geometry(), x );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TRANSFORMED_HH
