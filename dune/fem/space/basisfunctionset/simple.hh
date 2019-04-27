#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_SIMPLE_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_SIMPLE_HH

#include <cassert>
#include <cstddef>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/basisfunctionset/functor.hh>

namespace Dune
{

  namespace Fem
  {

    // SimpleBasisFunctionSet
    // ----------------------

    /** \class SimpleBasisFunctionSet
     *
     *  \brief This class is a simple basis function set which is needed for
     *  global basis functions sets (Fourier space etc.).
     *
     *  \note For localized basis function sets use the DefaultBasisFunctionSet.
     *
     *  \tparam  LocalFunctionSet  set of basis functions
     */
    template< class LocalFunctionSet >
    class SimpleBasisFunctionSet
    {
      typedef SimpleBasisFunctionSet< LocalFunctionSet > ThisType;

    public:
      typedef LocalFunctionSet LocalFunctionSetType;

      //! \brief entity type
      typedef typename LocalFunctionSetType::EntityType EntityType;
      typedef typename LocalFunctionSetType::FunctionSpaceType FunctionSpaceType;

      //! \brief range type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename FunctionSpaceType::RangeFieldType  RangeFieldType;
      typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

      //! \brief type of reference element
      typedef Dune::ReferenceElement< typename EntityType::Geometry > ReferenceElementType;

      /* default constructor
       *
       * Note: we require a local function set to have a default constructor;
       *       eventually use LocalFunctionSetProxy.
       */
      SimpleBasisFunctionSet () {}

      /** \brief constructor
       *
       *  \param[in]  localFunctionSet  implementation of LocalFunctionSet
       */
      explicit SimpleBasisFunctionSet ( const LocalFunctionSetType &localFunctionSet )
      : localFunctionSet_( localFunctionSet )
      {}

      //! \brief return order of basis function set
      int order () const { return localFunctionSet().order(); }

      //! \brief return size of basis function set
      std::size_t size () const { return localFunctionSet().size(); }

      //! \brief return reference element
      decltype(auto) referenceElement () const
      {
        return Dune::referenceElement< typename EntityType::Geometry::ctype, EntityType::Geometry::coorddimension >( entity().type() );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Quadrature, class Vector, class DofVector >
      void axpy ( const Quadrature &quad, const Vector &values, DofVector &dofs ) const
      {
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
      template< class Quadrature, class VectorA, class VectorB, class DofVector >
      void axpy ( const Quadrature &quad, const VectorA &valuesA, const VectorB &valuesB, DofVector &dofs ) const
      {
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
        FunctionalAxpyFunctor< RangeType, DofVector > functor( valueFactor, dofs );
        localFunctionSet().evaluateEach( x, functor );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        FunctionalAxpyFunctor< JacobianRangeType, DofVector > functor( jacobianFactor, dofs );
        localFunctionSet().jacobianEach( x, functor );
      }
      /** \brief Add H:D^2phi to each dof
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const HessianRangeType &hessianFactor, DofVector &dofs ) const
      {
        FunctionalAxpyFunctor< HessianRangeType, DofVector > functor( hessianFactor, dofs );
        localFunctionSet().hessianEach( x, functor );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor,
                  const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      /** \brief evaluate all basis functions and store the result in the
       *         ranges array
       */
      template< class Quadrature, class DofVector, class RangeArray >
      void evaluateAll ( const Quadrature &quad, const DofVector &dofs, RangeArray &ranges ) const
      {
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
      }

      //! please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        value = RangeType( 0 );
        AxpyFunctor< DofVector, RangeType > functor( dofs, value );
        localFunctionSet().evaluateEach( x, functor );
      }

      //! please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        AssignFunctor< RangeArray > functor( values );
        localFunctionSet().evaluateEach( x, functor );
      }

      //! please doc me
      template< class Quadrature, class DofVector, class JacobianRangeArray >
      void jacobianAll ( const Quadrature &quad, const DofVector &dofs, JacobianRangeArray &jacobians ) const
      {
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
      }

      //! please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        jacobian = JacobianRangeType( 0 );
        AxpyFunctor< DofVector, JacobianRangeType > functor( dofs, jacobian );
        localFunctionSet().jacobianEach( x, functor );
      }

      //! please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        AssignFunctor< JacobianRangeArray > functor( jacobians );
        localFunctionSet().jacobianEach( x, functor );
      }

      //! please doc me
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        hessian = HessianRangeType( typename HessianRangeType::value_type( typename RangeType::value_type( 0 ) ) );
        AxpyFunctor< DofVector, HessianRangeType > functor( dofs, hessian );
        localFunctionSet().hessianEach( x, functor );
      }

      //! please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        AssignFunctor< HessianRangeArray > functor( hessians );
        localFunctionSet().hessianEach( x, functor );
      }

      //! please doc me
      const EntityType &entity () const { return localFunctionSet().entity(); }


      // Non-interface methods
      // ---------------------

      //! \brief return local function set
      const LocalFunctionSetType localFunctionSet () const { return localFunctionSet_; }

    private:
      LocalFunctionSetType localFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_SIMPLE_HH
