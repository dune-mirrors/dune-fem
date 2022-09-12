#ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
#define DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH

// C++ includes
#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/shapefunctionset/caching.hh>

// dune-fem includes
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/basisfunctionset/transformation.hh>
#include <dune/fem/space/shapefunctionset/caching.hh>
#include <dune/fem/quadrature/cachingpointlist.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

#include <dune/fem/space/basisfunctionset/codegen.hh>
#include <dune/fem/space/basisfunctionset/evaluatecaller.hh>

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
     * \note ShapeFunctionSet must be a copyable object. For most
     *       non-trivial implementations, you may want to use a
     *       proxy, see file
\code
    <dune/fem/space/shapefunctionset/proxy.hh>
\endcode
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

      // if underlying shape function set was created with storage CodegenStorage
      // then this value should be true (see selectcaching.hh)
      static constexpr bool codegenShapeFunctionSet = detail::IsCodegenShapeFunctionSet< ShapeFunctionSetType >::value;

    protected:
      typedef typename ShapeFunctionSetType::FunctionSpaceType   LocalFunctionSpaceType;
      typedef typename LocalFunctionSpaceType::JacobianRangeType LocalJacobianRangeType;
      typedef typename LocalFunctionSpaceType::HessianRangeType  LocalHessianRangeType;

      typedef typename LocalFunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename EntityType::Geometry Geometry ;

      typedef typename Geometry::ctype ctype;
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
        : entity_( &entity ),
          shapeFunctionSet_( shapeFunctionSet )
      {
        // Note that this should be geometry_ = entity.geometry()
        // But Dune::Geometries are not assignable ...
        geometry_.reset();
        geometry_.emplace( entity.geometry() );
      }

      DefaultBasisFunctionSet ( const DefaultBasisFunctionSet &other )
        : entity_( other.entity_ ),
          shapeFunctionSet_( other.shapeFunctionSet_ )
      {
        // Note that this should be geometry_ = entity.geometry()
        // But Dune::Geometries are not assignable ...
        geometry_.reset();
        if( other.geometry_ )
          geometry_.emplace( other.geometry_.value() );
      }

      DefaultBasisFunctionSet &operator= ( const DefaultBasisFunctionSet &other )
      {
        entity_ = other.entity_;
        shapeFunctionSet_ = other.shapeFunctionSet_;

        // Note that this should be geometry_ = entity.geometry()
        // But Dune::Geometries are not assignable ...
        geometry_.reset();
        if( other.geometry_ )
          geometry_.emplace( other.geometry_.value() );
        return *this;
      }

      // Basis Function Set Interface Methods
      // ------------------------------------

      //! \brief return order of basis function set
      int order () const { return shapeFunctionSet().order(); }

      //! \brief return size of basis function set
      std::size_t size () const { return shapeFunctionSet().size(); }

      //! \brief return reference element
      auto referenceElement () const
        -> decltype( Dune::ReferenceElements< ctype, Geometry::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) )
      {
        return Dune::ReferenceElements< ctype, Geometry::coorddimension >::general( type() );
      }

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
        // don't know how to work directly with the DiagonalMatrix
        // returned from YaspGrid since gjit[k][l] is not correctly
        // implemented for k!=l so convert first into a dense matrix...
        Dune::FieldMatrix<double,DomainType::dimension,DomainType::dimension>
           G = gjit;
        DomainType Hg;
        for( int r = 0; r < FunctionSpaceType::dimRange; ++r )
          for( int j = 0; j < gjit.cols; ++j )
            for( int k = 0; k < gjit.cols; ++k )
            {
              hessianFactor[r].mv(G[k],Hg);
              tmpHessianFactor[r][j][k] += Hg * G[j];
            }
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
            const Geometry &geo = geometry();
            baseEval->evaluateJacobians( quad, geo, dofs, jacobians );
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
        const Geometry &geo = geometry();
        Transformation transformation( geo, coordinate( x ) );
        transformation( localJacobian, jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        assert( jacobians.size() >= size() );
        typedef JacobianTransformation< Geometry > Transformation;
        const Geometry &geo = geometry();

        Transformation transformation( geo, coordinate( x ) );
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
        const Geometry &geo = geometry();
        Transformation transformation( geo, coordinate( x ) );
        transformation( localHessian, hessian );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        assert( hessians.size() >= size() );
        typedef HessianTransformation< Geometry > Transformation;
        const Geometry &geo = geometry();
        Transformation transformation( geo, coordinate( x ) );
        AssignFunctor< HessianRangeArray, Transformation > f( hessians, transformation );
        shapeFunctionSet().hessianEach( x, f );
      }

      //! \brief return entity
      const Entity &entity () const
      {
        assert( valid() );
        return *entity_;
      }

      //! \brief return true if entity pointer is set
      bool valid () const { return bool(entity_); }

      //! \brief return geometry type
      Dune::GeometryType type () const { return entity().type(); }

      // Non-interface methods
      // ---------------------

      //! \brief return shape function set
      const ShapeFunctionSetType &shapeFunctionSet () const { return shapeFunctionSet_; }

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
            const Geometry &geo = geometry();
            baseEval->axpyJacobians( quad, geo, jacobianFactors, dofs );
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

    protected:
      Geometry geometry () const { return geometry_.value(); }

      const EntityType *entity_ = nullptr;
      ShapeFunctionSetType shapeFunctionSet_;
      std::optional< Geometry > geometry_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
