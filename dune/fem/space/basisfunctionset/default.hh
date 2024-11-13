#ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
#define DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH

// C++ includes
#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>
#include <optional>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/storage/entitygeometry.hh>
#include <dune/fem/space/shapefunctionset/caching.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/basisfunctionset/transformation.hh>
#include <dune/fem/space/shapefunctionset/caching.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

#include <dune/fem/space/basisfunctionset/codegen.hh>
#include <dune/fem/space/basisfunctionset/evaluatecaller.hh>

namespace Dune
{

  namespace Fem
  {
    template< class Entity, class ShapeFunctionSet >
    class BasisFunctionSetStorage : public EntityGeometryStorage< Entity >
    {
      typedef EntityGeometryStorage< Entity >  BaseType;
      typedef BasisFunctionSetStorage< Entity, ShapeFunctionSet > ThisType;

    public:
      //! \brief entity type
      typedef typename BaseType :: EntityType  EntityType;
      //! \brief type of geometry
      typedef typename BaseType :: Geometry    Geometry;

      //! \brief shape function set type
      typedef ShapeFunctionSet  ShapeFunctionSetType;

      // if underlying shape function set was created with storage CodegenStorage
      // then this value should be true (see selectcaching.hh)
      static constexpr bool codegenShapeFunctionSet = detail::IsCodegenShapeFunctionSet< ShapeFunctionSetType >::value;

      static const int pointSetId = detail::SelectPointSetId< ShapeFunctionSetType >::value;

      //! \brief constructor
      BasisFunctionSetStorage () {}

      //! \brief constructor
      explicit BasisFunctionSetStorage( const EntityType &entity, const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
        : BaseType( entity ),
          shapeFunctionSet_( shapeFunctionSet )
      {
      }

      BasisFunctionSetStorage ( const BasisFunctionSetStorage &other )
        : BaseType( other ),
          shapeFunctionSet_( other.shapeFunctionSet_ )
      {
      }

      BasisFunctionSetStorage &operator= ( const BasisFunctionSetStorage &other )
      {
        BaseType::operator=(other);
        shapeFunctionSet_ = other.shapeFunctionSet_;
        return *this;
      }

      using BaseType :: entity;
      using BaseType :: valid;
      using BaseType :: geometry;
      using BaseType :: type;
      using BaseType :: referenceElement;

      // Non-interface methods
      // ---------------------

      //! \brief return shape function set
      const ShapeFunctionSetType &shapeFunctionSet () const { return shapeFunctionSet_; }

      // Basis Function Set Interface Methods
      // ------------------------------------

      //! \brief return order of basis function set
      int order () const { return shapeFunctionSet().order(); }

      //! \brief return size of basis function set
      std::size_t size () const { return shapeFunctionSet().size(); }

    protected:
      ShapeFunctionSetType shapeFunctionSet_;
    };


    // DefaultBasisFunctionSet
    // -----------------------

    /**
     * \brief implementation of a basis function set for given entity
     *
     * \tparam  Entity            entity type
     * \tparam  ShapeFunctionSet  shape function set
     *
     * \note ShapeFunctionSet must be a copyable object. For most
     *       non-trivial implementations, you may want to use a
     *       proxy, see file
\code
    <dune/fem/space/shapefunctionset/proxy.hh>
\endcode
     */
    template< class Entity, class ShapeFunctionSet >
    class DefaultBasisFunctionSet
      : public BasisFunctionSetStorage< Entity, ShapeFunctionSet >
    {
      typedef BasisFunctionSetStorage< Entity, ShapeFunctionSet > BaseType;
      typedef DefaultBasisFunctionSet< Entity, ShapeFunctionSet > ThisType;

    public:
      //! \brief entity type
      typedef typename BaseType::EntityType  EntityType;

      //! \brief geometry
      typedef typename BaseType::Geometry    Geometry ;

      //! \brief type of coordinate field
      typedef typename Geometry::ctype       ctype;

      //! \brief shape function set type
      typedef typename BaseType::ShapeFunctionSetType  ShapeFunctionSetType;

      // if underlying shape function set was created with storage CodegenStorage
      // then this value should be true (see selectcaching.hh)
      using BaseType :: codegenShapeFunctionSet;

    protected:
      typedef typename ShapeFunctionSetType::FunctionSpaceType   LocalFunctionSpaceType;
      typedef typename LocalFunctionSpaceType::JacobianRangeType LocalJacobianRangeType;
      typedef typename LocalFunctionSpaceType::HessianRangeType  LocalHessianRangeType;

      typedef typename LocalFunctionSpaceType::RangeFieldType RangeFieldType;

    public:
      //  slight misuse of struct ToLocalFunctionSpace!!!
      //! \brief type of function space
      typedef typename ToNewDimDomainFunctionSpace< LocalFunctionSpaceType, Geometry::coorddimension > :: Type  FunctionSpaceType;

      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      typedef typename FunctionSpaceType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

      typedef typename ScalarFunctionSpaceType::RangeType         ScalarRangeType;
      typedef typename ScalarFunctionSpaceType::JacobianRangeType ScalarJacobianRangeType;

      //! \brief type of reference element
      typedef std::decay_t< decltype( Dune::ReferenceElements< ctype, Geometry::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) ) > ReferenceElementType;

      static const int pointSetId = detail::SelectPointSetId< ShapeFunctionSetType >::value;

      //! \brief constructor
      DefaultBasisFunctionSet () {}

      //! \brief constructor
      explicit DefaultBasisFunctionSet ( const EntityType &entity, const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
        : BaseType( entity, shapeFunctionSet )
      {
      }

      DefaultBasisFunctionSet ( const DefaultBasisFunctionSet &other ) = default;
      DefaultBasisFunctionSet &operator= ( const DefaultBasisFunctionSet &other ) = default;

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
        axpyImpl( quad, values, dofs, values[ 0 ] );
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
        assert( valuesA.size() > 0 );
        assert( valuesB.size() > 0 );

        axpyImpl( quad, valuesA, dofs, valuesA[ 0 ] );
        axpyImpl( quad, valuesB, dofs, valuesB[ 0 ] );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        FunctionalAxpyFunctor< RangeType, DofVector > f( valueFactor, dofs );
        shapeFunctionSet().evaluateEach( x, f );
      }

      /** \brief evaluate all derivatives of all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        typedef typename Geometry::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const Geometry &geo = geometry();
        const GeometryJacobianInverseTransposedType &gjit = geo.jacobianInverseTransposed( coordinate( x ) );
        LocalJacobianRangeType tmpJacobianFactor;
        for( int r = 0; r < FunctionSpaceType::dimRange; ++r )
          gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

        FunctionalAxpyFunctor< LocalJacobianRangeType, DofVector > f( tmpJacobianFactor, dofs );
        shapeFunctionSet().jacobianEach( x, f );
      }

      /** \brief evaluate all basis function and derivatives and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      /** \brief Add H:D^2phi to each dof
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const HessianRangeType &hessianFactor, DofVector &dofs ) const
      {
        typedef typename Geometry::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const Geometry &geo = geometry();
        const GeometryJacobianInverseTransposedType &gjit = geo.jacobianInverseTransposed( coordinate( x ) );
        LocalHessianRangeType tmpHessianFactor( RangeFieldType(0) );
        hessianTransformation(gjit.transposed(),hessianFactor,tmpHessianFactor);

        FunctionalAxpyFunctor< LocalHessianRangeType, DofVector > f( tmpHessianFactor, dofs );
        shapeFunctionSet().hessianEach( x, f );
      }


      /** \copydoc BasisFunctionSet::evaluateAll( quad, dofs, ranges ) */
      template< class QuadratureType, class DofVector, class RangeArray >
      void evaluateAll ( const QuadratureType &quad, const DofVector &dofs, RangeArray &ranges ) const
      {
        assert( ranges.size() >= quad.nop() );

        // if shape function set supports codegen and quadrature supports caching
        if constexpr ( codegenShapeFunctionSet && std::is_base_of< CachingInterface, QuadratureType > :: value)
        {
          typedef Codegen :: EvaluateCallerInterfaceTraits< QuadratureType, RangeArray, DofVector > Traits;
          typedef Codegen :: EvaluateCallerInterface< Traits > BaseEvaluationType;

          // get base function evaluate caller
          const auto& baseEval = BaseEvaluationType::storage( *this, rangeCache( quad ), quad );

          // true if implementation exists, otherwise this is a nullptr
          if( baseEval )
          {
            baseEval->evaluateRanges( quad, dofs, ranges );
            return ;
          }
        }

        {
          // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
          const unsigned int nop = quad.nop();
          for( unsigned int qp = 0; qp < nop; ++qp )
          {
            evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
          }
        }
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        value = RangeType( 0 );
        AxpyFunctor< DofVector, RangeType > f( dofs, value );
        shapeFunctionSet().evaluateEach( x, f );
      }

      //! \todo please doc me
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        assert( values.size() >= size() );
        std::fill( values.begin(), values.end(), RangeType(0) );
        AssignFunctor< RangeArray > f( values );
        shapeFunctionSet().evaluateEach( x, f );
      }

      /** \copydoc BasisFunctionSet::jacobianAll( quad, dofs, jacobians ) */
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        assert( jacobians.size() >= quad.nop() );

        // if shape function set supports codegen and quadrature supports caching
        if constexpr ( codegenShapeFunctionSet && std::is_base_of< CachingInterface, QuadratureType > :: value)
        {
          typedef Codegen :: EvaluateCallerInterfaceTraits< QuadratureType, JacobianArray, DofVector, Geometry >  Traits;
          typedef Codegen :: EvaluateCallerInterface< Traits > BaseEvaluationType;

          // get base function evaluate caller (calls axpyRanges)
          const auto& baseEval = BaseEvaluationType::storage( *this, jacobianCache( quad ), quad );

          // true if implementation exists
          if( baseEval )
          {
            // call appropriate axpyJacobian method
            baseEval->evaluateJacobians( quad, geometry(), dofs, jacobians );
            return ;
          }
        }

        {
          // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
          const unsigned int nop = quad.nop();
          for( unsigned int qp = 0; qp < nop; ++qp )
          {
            jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
          }
        }
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        LocalJacobianRangeType localJacobian( RangeFieldType( 0 ) );
        AxpyFunctor< DofVector, LocalJacobianRangeType > f( dofs, localJacobian );
        shapeFunctionSet().jacobianEach( x, f );

        typedef JacobianTransformation< Geometry > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        transformation( localJacobian, jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        assert( jacobians.size() >= size() );
        typedef JacobianTransformation< Geometry > Transformation;

        Transformation transformation( geometry(), coordinate( x ) );
        AssignFunctor< JacobianRangeArray, Transformation > f( jacobians, transformation );
        shapeFunctionSet().jacobianEach( x, f );
      }

      //! \todo please doc me
      template< class QuadratureType, class DofVector, class HessianArray >
      void hessianAll ( const QuadratureType &quad, const DofVector &dofs, HessianArray &hessians ) const
      {
        assert( hessians.size() >= quad.nop() );
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          hessianAll( quad[ qp ], dofs, hessians[ qp ] );
        }
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        LocalHessianRangeType localHessian( typename LocalHessianRangeType::value_type( RangeFieldType( 0 ) ) );
        AxpyFunctor< DofVector, LocalHessianRangeType > f( dofs, localHessian );
        shapeFunctionSet().hessianEach( x, f );

        typedef HessianTransformation< Geometry > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        transformation( localHessian, hessian );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        assert( hessians.size() >= size() );
        typedef HessianTransformation< Geometry > Transformation;
        Transformation transformation( geometry(), coordinate( x ) );
        AssignFunctor< HessianRangeArray, Transformation > f( hessians, transformation );
        shapeFunctionSet().hessianEach( x, f );
      }
    protected:
      //! \brief evaluate all basis function and multiply with given values and add to dofs
      template< class QuadratureType, class RangeArray, class DofVector >
      void axpyImpl ( const QuadratureType &quad, const RangeArray &rangeFactors, DofVector &dofs, const RangeType& ) const
      {
        assert( rangeFactors.size() >= quad.nop() );

        // if shape function set supports codegen and quadrature supports caching
        if constexpr ( codegenShapeFunctionSet && std::is_base_of< CachingInterface, QuadratureType > :: value)
        {
          typedef Codegen :: EvaluateCallerInterfaceTraits< QuadratureType, RangeArray, DofVector > Traits;
          typedef Codegen :: EvaluateCallerInterface< Traits > BaseEvaluationType;

          // get base function evaluate caller (calls axpyRanges)
          const auto& baseEval = BaseEvaluationType::storage( *this, rangeCache( quad ), quad );

          // true if implementation exists
          if( baseEval )
          {
            baseEval->axpyRanges( quad, rangeFactors, dofs );
            return ;
          }
        }

        {
          // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
          const unsigned int nop = quad.nop();
          for( unsigned int qp = 0; qp < nop; ++qp )
          {
            axpy( quad[ qp ], rangeFactors[ qp ], dofs );
          }
        }
      }

      //! \brief evaluate all basis function and multiply with given values and add to dofs
      template< class QuadratureType, class JacobianArray, class DofVector >
      void axpyImpl ( const QuadratureType &quad, const JacobianArray &jacobianFactors, DofVector &dofs, const JacobianRangeType& ) const
      {
        assert( jacobianFactors.size() >= quad.nop() );
        // if shape function set supports codegen and quadrature supports caching
        if constexpr ( codegenShapeFunctionSet && std::is_base_of< CachingInterface, QuadratureType > :: value)
        {
          typedef Codegen :: EvaluateCallerInterfaceTraits< QuadratureType, JacobianArray, DofVector, Geometry >  Traits;
          typedef Codegen :: EvaluateCallerInterface< Traits > BaseEvaluationType;

          // get base function evaluate caller (calls axpyRanges)
          const auto& baseEval = BaseEvaluationType::storage( *this, jacobianCache( quad ), quad );

          // true if implementation exists
          if( baseEval )
          {
            // call appropriate axpyRanges method
            baseEval->axpyJacobians( quad, geometry(), jacobianFactors, dofs );
            return ;
          }
        }

        {
          // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
          const unsigned int nop = quad.nop();
          for( unsigned int qp = 0; qp < nop; ++qp )
          {
            axpy( quad[ qp ], jacobianFactors[ qp ], dofs );
          }
        }
      }

      //! \brief evaluate all basis function and multiply with given values and add to dofs
      template< class QuadratureType, class HessianArray, class DofVector >
      void axpyImpl ( const QuadratureType &quad, const HessianArray &hessianFactors, DofVector &dofs, const HessianRangeType& ) const
      {
        assert( hessianFactors.size() >= quad.nop() );
        /* TODO add code generation for hessians
        // if shape function set supports codegen and quadrature supports caching
        if constexpr ( codegenShapeFunctionSet && std::is_base_of< CachingInterface, QuadratureType > :: value)
        {
          typedef Codegen :: EvaluateCallerInterfaceTraits< QuadratureType, HessianArray, DofVector, Geometry >  Traits;
          typedef Codegen :: EvaluateCallerInterface< Traits > BaseEvaluationType;

          // get base function evaluate caller (calls axpyRanges)
          const auto& baseEval = BaseEvaluationType::storage( *this, hessianCache( quad ), quad );

          // true if implementation exists
          if( baseEval )
          {
            // call appropriate axpyRanges method
            const Geometry &geo = geometry();
            baseEval->axpyHessian( quad, geo, hessianFactors, dofs );
            return ;
          }
        }
        */
        {
          const unsigned int nop = quad.nop();
          for( unsigned int qp = 0; qp < nop; ++qp )
          {
            axpy( quad[ qp ], hessianFactors[ qp ], dofs );
          }
        }
      }

      template <class QuadratureType>
      const auto& rangeCache( const QuadratureType& quad ) const
      {
        return shapeFunctionSet().scalarShapeFunctionSet().impl().rangeCache( quad );
      }

      template <class QuadratureType>
      const auto& jacobianCache( const QuadratureType& quad ) const
      {
        return shapeFunctionSet().scalarShapeFunctionSet().impl().jacobianCache( quad );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
