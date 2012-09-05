#ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
#define DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH

//- C++ includes
#include <cassert>
#include <cstddef>

//- dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

//- dune-fem includes
#include <dune/fem/space/basefunctionset/functor.hh>
#include <dune/fem/space/basefunctionset/transformation.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    // DefaultBasisFunctionSet
    // -----------------------

    /**
     * \brief implementation of a basis function set for given entity
     *
     * \tparam  Entity            entity type
     * \tparam  ShapeFunctionSet  shape function set
     *
     */
    template< class Entity, class ShapeFunctionSet >
    class DefaultBasisFunctionSet
    {
      typedef DefaultBasisFunctionSet< Entity, ShapeFunctionSet > ThisType;

    public:
      //! \brief entity type
      typedef Entity EntityType;
      //! \brief shape function set type
      typedef ShapeFunctionSet ShapeFunctionSetType;

    protected:
      typedef typename ShapeFunctionSetType::FunctionSpaceType LocalFunctionSpaceType;
      typedef typename LocalFunctionSpaceType::JacobianRangeType LocalJacobianRangeType;
      typedef typename LocalFunctionSpaceType::HessianRangeType LocalHessianRangeType;

      typedef typename LocalFunctionSpaceType::RangeFieldType RangeFieldType;

    public:
      //  slight misuse of struct ToLocalFunctionSpace!!!
      //! \brief type of function space
      typedef typename ToLocalFunctionSpace< LocalFunctionSpaceType, EntityType::Geometry::coorddimension >::Type FunctionSpaceType;

      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief type of reference element
      typedef Dune::ReferenceElement< typename EntityType::Geometry::ctype, 
                                      EntityType::Geometry::coorddimension > ReferenceElementType;

      //! \brief constructor
      DefaultBasisFunctionSet ( const EntityType &entity, const ShapeFunctionSet &shapeFunctionSet )
      : entity_( &entity ),
        shapeFunctionSet_( shapeFunctionSet )
      {
        assert( entity.type() == shapeFunctionSet.type() );
      }


      // Basis Function Set Interface Methods
      // ------------------------------------

      //! \brief return size of basis function set
      std::size_t size () const { return shapeFunctionSet().size(); }

      //! \brief return geometry type of entity
      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const { return shapeFunctionSet().type(); }

      //! \brief return reference element
      const ReferenceElementType &referenceElement () const
      {
        return Dune::ReferenceElements< typename EntityType::Geometry::ctype, 
                                        EntityType::Geometry::coorddimension >::general( shapeFunctionSet().type() );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        FunctionalAxpyFunctor< RangeType, DofVector > f( valueFactor, dofs );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        typedef typename GeometryType::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const GeometryJacobianInverseTransposedType &gjit = geometry().jacobianInverseTransposed( coordinate( x ) );
        LocalJacobianRangeType tmpJacobianFactor;
        for( int r = 0; r < FunctionSpaceType::dimRange; ++r )
          gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

        FunctionalAxpyFunctor< LocalJacobianRangeType, DofVector > f( jacobianFactor, dofs );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        AxpyFunctor< DofVector, RangeType > f( dofs, value );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        AssignFunctor< RangeArray > f( values );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        LocalJacobianRangeType localJacobian( RangeFieldType( 0 ) );
        AxpyFunctor< DofVector, LocalJacobianRangeType > f( dofs, localJacobian );
        shapeFunctionSet().jacobianEach( x, f );

        typedef JacobianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        transformation( localJacobian, jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        typedef JacobianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        AssignFunctor< JacobianRangeArray, Transformation > f( jacobians, transformation );
        shapeFunctionSet().jacobianEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        LocalHessianRangeType localHessian( RangeFieldType( 0 ) );
        AxpyFunctor< DofVector, LocalHessianRangeType > f( dofs, localHessian );
        shapeFunctionSet().hessianEach( x, f );

        typedef HessianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        transformation( localHessian, hessian );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        typedef HessianTransformation< GeometryType > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        AssignFunctor< HessianRangeArray, Transformation > f( hessians, transformation );
        shapeFunctionSet().hessianEach( x, f );
      }

      //! \brief return entity
      const Entity &entity () const { return *entity_; }


      // Non-interface methods
      // ---------------------

      //! \brief return shape function set
      const ShapeFunctionSet &shapeFunctionSet () const { return shapeFunctionSet_; }

    protected:
      typedef typename EntityType::Geometry GeometryType;
      const GeometryType &geometry () const { return entity().geometry(); }

    private:
      const EntityType *entity_;
      const ShapeFunctionSetType &shapeFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
