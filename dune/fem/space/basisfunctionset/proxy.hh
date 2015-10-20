#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_PROXY_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_PROXY_HH

// C++ includes
#include <cassert>
#include <cstddef>

// dune-geometry includes
#include <dune/geometry/type.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Provides a proxy class for pointers to a basis function set
*/


namespace Dune
{

  namespace Fem
  {

    // BasisFunctionSetProxy
    // ---------------------

    /*
     * \brief A proxy object converting a pointer to a basis function set to a object
     *
     * \tparam  BasisFunctionSet  An implementation of Dune::Fem::BasisFunctionSet
     *
     * \note This class has an implicit constructor from a pointer to a basis function set.
     */

    template< class BasisFunctionSet >
    class BasisFunctionSetProxy
    {
      typedef BasisFunctionSetProxy< BasisFunctionSet > ThisType;

    public:
      typedef BasisFunctionSet ImplementationType;
      const ImplementationType &impl () const
      {
        assert( basisFunctionSet_ );
        return *basisFunctionSet_;
      }

      typedef typename BasisFunctionSet::EntityType EntityType;

      typedef typename BasisFunctionSet::FunctionSpaceType FunctionSpaceType;

      typedef typename BasisFunctionSet::DomainType DomainType;
      typedef typename BasisFunctionSet::RangeType RangeType;
      typedef typename BasisFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename BasisFunctionSet::HessianRangeType HessianRangeType;

      typedef typename BasisFunctionSet::ReferenceElementType ReferenceElementType;

      BasisFunctionSetProxy ()
      : basisFunctionSet_( nullptr )
      {}

      BasisFunctionSetProxy ( const BasisFunctionSet *basisFunctionSet )
      : basisFunctionSet_( basisFunctionSet )
      {}

      //! \copydoc BasisFunctionSet::order
      int order () const { return impl().order(); }

      //! \copydoc BasisFunctionSet::size
      std::size_t size () const { return impl().size(); }

      //! \copydoc BasisFunctionSet::type
      Dune::GeometryType type () const { return impl().type(); }

      //! \copydoc BasisFunctionSet::referenceElement
      const ReferenceElementType &referenceElement () const
      {
        return impl().referenceElement();
      }

      //! \copydoc BasisFunctionSet::axpy( quad, values, dofs )
      template< class QuadratureType, class Vector, class DofVector >
      void axpy ( const QuadratureType &quad, const Vector &values, DofVector &dofs ) const
      {
        impl().axpy( quad, values, dofs );
      }

      //! \copydoc BasisFunctionSet::axpy( quad, valuesA, valuesB, dofs )
      template< class QuadratureType, class VectorA, class VectorB, class DofVector >
      void axpy ( const QuadratureType &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
      {
        impl().axpy( quad, valuesA, valuesB, dofs );
      }

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        impl().axpy( x, valueFactor, dofs );
      }

      //! \copydoc BasisFunctionSet::axpy( x, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        axpy( x, jacobianFactor, dofs );
      }

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        impl().axpy( x, valueFactor, jacobianFactor, dofs );
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, dofs, ranges )
      template< class Point, class DofVector, class RangeArray  >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeArray& ranges ) const
      {
        impl().evaluateAll( x, dofs, ranges );
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, dofs, value )
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        impl().evaluateAll( x, dofs, value );
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, values )
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        impl().evaluateAll( x, values );
      }

      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobians )
      template< class Point, class DofVector, class JacobianArray >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        impl().jacobianAll( x, dofs, jacobians );
      }

      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobian )
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        impl().jacobianAll( x, dofs, jacobian );
      }

      //! \copydoc BasisFunctionSet::jacobianAll( x, jacobians )
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        impl().jacobianAll( x, jacobians );
      }

      //! \copydoc BasisFunctionSet::hessianAll( x, dofs, hessian )
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        impl().hessianAll( x, dofs, hessian );
      }

      //! \copydoc BasisFunctionSet::hessianAll( x, hessians )
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        impl().hessianAll( x, hessians );
      }

      //! \copydoc BasisFunctionSet::entity
      const EntityType &entity () const { return impl().entity(); }

    private:
      const BasisFunctionSet *basisFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_PROXY_HH
