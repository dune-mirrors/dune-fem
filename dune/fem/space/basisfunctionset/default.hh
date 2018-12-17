#ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
#define DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH


// C++ includes
#include <cassert>
#include <cstddef>
#include <memory>

// dune-common includes
#include <dune/common/std/optional.hh>

#include <type_traits>
#include <utility>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/common/coordinate.hh>
#include <dune/fem/space/basisfunctionset/functor.hh>
#include <dune/fem/space/basisfunctionset/transformation.hh>
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

    protected:
      typedef typename ShapeFunctionSetType::FunctionSpaceType LocalFunctionSpaceType;
      typedef typename LocalFunctionSpaceType::JacobianRangeType LocalJacobianRangeType;
      typedef typename LocalFunctionSpaceType::HessianRangeType LocalHessianRangeType;

      typedef typename LocalFunctionSpaceType::RangeFieldType RangeFieldType;

      typedef typename EntityType::Geometry GeometryType;

      typedef typename GeometryType::ctype ctype;

    public:
      //  slight misuse of struct ToLocalFunctionSpace!!!
      //! \brief type of function space
      typedef typename ToNewDimDomainFunctionSpace< LocalFunctionSpaceType, EntityType :: Geometry :: coorddimension >::Type FunctionSpaceType;

      //! \brief domain type
      typedef typename FunctionSpaceType::DomainType DomainType;

      //! \brief range type
      typedef typename FunctionSpaceType::RangeType RangeType;
      //! \brief jacobian range type
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      //! \brief hessian range type
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      //! \brief type of reference element
      typedef std::decay_t< decltype( Dune::ReferenceElements< ctype, GeometryType::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) ) > ReferenceElementType;

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
        -> decltype( Dune::ReferenceElements< ctype, GeometryType::coorddimension >::general( std::declval< const Dune::GeometryType & >() ) )
      {
        return Dune::ReferenceElements< ctype, GeometryType::coorddimension >::general( type() );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class QuadratureType, class Vector, class DofVector >
      void axpy ( const QuadratureType &quad, const Vector &values, DofVector &dofs ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          axpy( quad[ qp ], values[ qp ], dofs );
        }
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
        FunctionalAxpyFunctor< RangeType, DofVector > f( valueFactor, dofs );
        shapeFunctionSet().evaluateEach( x, f );
      }

      /** \brief evaluate all basis function and multiply with given
       *         values and add to dofs
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        typedef typename GeometryType::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const GeometryType &geo = geometry();
        const GeometryJacobianInverseTransposedType &gjit = geo.jacobianInverseTransposed( coordinate( x ) );
        LocalJacobianRangeType tmpJacobianFactor( RangeFieldType( 0 ) );
        for( int r = 0; r < FunctionSpaceType::dimRange; ++r )
          gjit.mtv( jacobianFactor[ r ], tmpJacobianFactor[ r ] );

        FunctionalAxpyFunctor< LocalJacobianRangeType, DofVector > f( tmpJacobianFactor, dofs );
        shapeFunctionSet().jacobianEach( x, f );
      }
      /** \brief Add H:D^2phi to each dof
       */
      template< class Point, class DofVector >
      void axpy ( const Point &x, const HessianRangeType &hessianFactor, DofVector &dofs ) const
      {
        typedef typename GeometryType::JacobianInverseTransposed GeometryJacobianInverseTransposedType;
        const GeometryType &geo = geometry();
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
        {
          evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
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
        AssignFunctor< RangeArray > f( values );
        shapeFunctionSet().evaluateEach( x, f );
      }

      /** \copydoc BasisFunctionSet::jacobianAll( quad, dofs, jacobians ) */
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
        {
          jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
        }
      }

      //! \todo please doc me
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        LocalJacobianRangeType localJacobian( RangeFieldType( 0 ) );
        AxpyFunctor< DofVector, LocalJacobianRangeType > f( dofs, localJacobian );
        shapeFunctionSet().jacobianEach( x, f );
        const GeometryType &geo = geometry();
        typedef JacobianTransformation< GeometryType > Transformation;
        Transformation transformation( geo, coordinate( x ) );
        transformation( localJacobian, jacobian );
      }

      //! \todo please doc me
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        assert( jacobians.size() >= size() );
        const GeometryType &geo = geometry();
        typedef JacobianTransformation< GeometryType > Transformation;
        Transformation transformation( geo, coordinate( x ) );
        AssignFunctor< JacobianRangeArray, Transformation > f( jacobians, transformation );
        shapeFunctionSet().jacobianEach( x, f );
      }

      //! \todo please doc me
      template< class QuadratureType, class DofVector, class HessianArray >
      void hessianAll ( const QuadratureType &quad, const DofVector &dofs, HessianArray &hessians ) const
      {
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
        const GeometryType &geo = geometry();
        typedef HessianTransformation< GeometryType > Transformation;
        Transformation transformation( geo, coordinate( x ) );
        transformation( localHessian, hessian );
      }

      //! \todo please doc me
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        assert( hessians.size() >= size() );
        const GeometryType &geo = geometry();
        typedef HessianTransformation< GeometryType > Transformation;
        Transformation transformation( geo, coordinate( x ) );
        AssignFunctor< HessianRangeArray, Transformation > f( hessians, transformation );
        shapeFunctionSet().hessianEach( x, f );
      }

      //! \brief return entity
      const Entity &entity () const
      {
        assert( entity_ );
        return *entity_;
      }

      //! \brief return geometry type
      Dune::GeometryType type () const { return entity().type(); }


      // Non-interface methods
      // ---------------------

      //! \brief return shape function set
      const ShapeFunctionSetType &shapeFunctionSet () const { return shapeFunctionSet_; }

    protected:
      GeometryType geometry () const { return geometry_.value(); }

    private:
      const EntityType *entity_ = nullptr;
      ShapeFunctionSetType shapeFunctionSet_;
      Std::optional< GeometryType > geometry_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_DEFAULT_HH
