#ifndef DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH
#define DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH

//- C++ includes
#include <cstddef>

//- dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

//- dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>


namespace Dune
{

  namespace Fem
  {

    // BasisFunctionSet
    // ----------------

    template< class Entity, class Range, class Implementation >
    class BasisFunctionSet
    {
    public:
      //! \brief entity type
      typedef Entity EntityType;

      //! \brief implementation type
      typedef Implementation ImplementationType;
      //! \brief return reference to implementation
      ImplementationType &impl () { return impl_; }
      //! \brief return reference to implementation
      const ImplementationType &impl () const { return impl_; }

      // function space type
      typedef FunctionSpace< typename Entity::Geometry::ctype, typename Range::value_type, 
                             Entity::Geometry::coorddimension, Range::dimension
                           > FunctionSpaceType;
       
      //! \brief range type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief type of reference element
      typedef Dune::ReferenceElement< typename DomainType::value_type, 
                                      DomainType::dimension > ReferenceElementType;

      //! \brief copy constructor from implementation
      BasisFunctionSet ( const ImplementationType &impl )
      : impl_( impl )
      {}

      //! \brief return size of basis function set
      std::size_t size () const { return impl().size(); }

      //! \brief return geometry type of entity
      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const { return referenceElement().type(); }

      //! \brief return reference element
      const ReferenceElementType &referenceElement () const { return impl().referenceElement(); }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        impl().axpy( x, valueFactor, dofs );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        axpy( x, jacobianFactor, dofs );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        impl().axpy( x, valueFactor, jacobianFactor, dofs );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        impl().evaluateAll( x, dofs, value );
      }

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        impl().evaluateAll( x, values );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        impl().jacobianAll( x, dofs, jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        impl().jacobianAll( x, jacobians );
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        impl().hessianAll( x, dofs, hessian );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        impl().hessianAll( x, hessians );
      }

      //! \brief return entity
      const EntityType &entity () const { return impl().entity(); }

    private:
      ImplementationType impl_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH
