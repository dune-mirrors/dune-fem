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

      //! \brief return size of basis function set
      std::size_t size () const;

      //! \brief return geometry type of entity
      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const;

      //! \brief return reference element
      const ReferenceElementType &referenceElement () const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const;

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const;

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const;

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const;

      //! \brief return entity
      const EntityType &entity () const;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH
