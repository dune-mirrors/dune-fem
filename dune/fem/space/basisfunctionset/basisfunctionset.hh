#ifndef DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH
#define DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH

#include <cstddef>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/common/functionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // BasisFunctionSet
    // ----------------

    /**
     * \brief Interface class for basis function sets
     *
     * This class cannot be used itself, it is for documentation purposes
     * only.
     *
     * \note Constructor signatures are explicitly not specified by this
     *       interface.
     */
    template< class Entity, class Range >
    class BasisFunctionSet
    {
    public:
      //! \brief entity type
      typedef Entity EntityType;

      //! function space type
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

      //! \brief return order of basis function set
      int order () const;

      //! \brief return size of basis function set
      std::size_t size () const;

      //! \brief return reference element
      const ReferenceElementType &referenceElement () const;

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Quadrature, class Vector, class DofVector >
      void axpy ( const Quadrature &quad, const Vector &values, DofVector &dofs ) const;

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Quadrature, class VectorA, class VectorB, class DofVector >
      void axpy ( const Quadrature &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const;

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const;

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const;

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const;

      /** \brief evaluate all basis functions and store the result in the
       *         ranges array
       */
      template< class Quadrature, class DofVector, class RangeArray >
      void evaluateAll ( const Quadrature &quad, const DofVector &dofs, RangeArray &ranges ) const;

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const;

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const;

      //! \brief evaluate the jacobian of all basis functions and store the result in the jacobians array
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const;

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

      //! \brief return true if entity was set
      bool valid () const;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_BASISFUNCTIONSET_HH
