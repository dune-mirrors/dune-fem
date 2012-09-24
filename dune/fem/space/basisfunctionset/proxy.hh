#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_PROXY_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_PROXY_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/nullptr.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/basisfunctionset.hh>

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
    : public Dune::Fem::BasisFunctionSet< typename BasisFunctionSet::EntityType, 
                                          typename BasisFunctionSet::RangeType, 
                                          BasisFunctionSetProxy< BasisFunctionSet > > 
    {
      typedef BasisFunctionSetProxy< BasisFunctionSet > ThisType;
      typedef Dune::Fem::BasisFunctionSet< typename BasisFunctionSet::EntityType,
                                           typename BasisFunctionSet::RangeType,
                                           BasisFunctionSetProxy< BasisFunctionSet > > BaseType;

    public:
      typedef typename BaseType::EntityType EntityType;

      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::HessianRangeType HessianRangeType;

      typedef typename BaseType::ReferenceElementType ReferenceElementType; 

      BasisFunctionSetProxy ()
      : basisFunctionSet_( nullptr )
      {}
      
      BasisFunctionSetProxy ( const BasisFunctionSet *basisFunctionSet )
      : basisFunctionSet_( basisFunctionSet )
      {}

      std::size_t size () const { return basisFunctionSet().size(); } 

      Dune::GeometryType type () const { return basisFunctionSet().type(); } 

      const ReferenceElementType &referenceElement () const
      {
        return basisFunctionSet().referenceElement();
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        basisFunctionSet().axpy( x, valueFactor, dofs );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        axpy( x, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        basisFunctionSet().axpy( x, valueFactor, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        basisFunctionSet().evaluateAll( x, dofs, value );
      }

      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        basisFunctionSet().evaluateAll( x, values );
      }

      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        basisFunctionSet().jacobianAll( x, dofs, jacobian );
      }

      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        basisFunctionSet().jacobianAll( x, jacobians );
      }

      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        basisFunctionSet().hessianAll( x, dofs, hessian );
      }

      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        basisFunctionSet().hessianAll( x, hessians );
      }

      const EntityType &entity () const { return basisFunctionSet().entity(); }

    private:
      const BasisFunctionSet &basisFunctionSet () const
      {
        assert( basisFunctionSet_ );
        return *basisFunctionSet_;
      }

      BasisFunctionSet basisFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_PROXY_HH
