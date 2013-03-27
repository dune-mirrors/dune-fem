#ifndef DUNE_FEM_BASISFUNCTIONSET_VECTORIAL_HH
#define DUNE_FEM_BASISFUNCTIONSET_VECTORIAL_HH

#include <cstddef>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/version.hh>

namespace Dune
{

  namespace Fem
  {

    // VectorialBasisFunctionSet
    // -------------------------

    template< class ScalarBasisFunctionSet, class Range >
    class VectorialBasisFunctionSet
    {
      typedef VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range > ThisType;

    public:
      typedef ScalarBasisFunctionSet ScalarBasisFunctionSetType;

      typedef typename ScalarBasisFunctionSetType::EntityType EntityType;;
      typedef typename ScalarBasisFunctionSetType::ReferenceElementType ReferenceElementType;

    private:
      typedef typename ScalarBasisFunctionSetType::FunctionSpaceType ScalarFunctionSpaceType;
      static const int dimRange = Range::dimension;

    public:
      typedef typename ToNewDimRangeFunctionSpace< ScalarFunctionSpaceType, dimRange >::Type FunctionSpaceType;
       
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    private:
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      struct Evaluate;
      struct Jacobian;
      struct Hessian;

    public:
      explicit VectorialBasisFunctionSet ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet )
      : scalarBasisFunctionSet_( scalarBasisFunctionSet )
      {}

      int order () const { return scalarBasisFunctionSet().order(); }

      std::size_t size () const { return dimRange*scalarBasisFunctionSet().size(); }

      DUNE_VERSION_DEPRECATED(1,4,remove)
      Dune::GeometryType type () const { return scalarBasisFunctionSet().type(); }

      const ReferenceElementType &referenceElement () const { return scalarBasisFunctionSet().referenceElement(); }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        axpy< Evaluate >( x, valueFactor, dofs );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        axpy< Jacobian >( x, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor,
                  DofVector &dofs ) const
      {
        axpy( x, valueFactor, dofs );
        axpy( x, jacobianFactor, dofs );
      }

      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        evaluateAll< Evaluate >( x, dofs, value );
      }

      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        evaluateAll< Evaluate >( x, values );
      }

      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        evaluateAll< Jacobian >( x, dofs, jacobian );
      }

      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        evaluateAll< Jacobian >( x, jacobians );
      }

      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        evaluateAll< Hessian >( x, dofs, hessian );
      }

      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        evaluateAll< Hessian >( x, hessians );
      }

      const EntityType &entity () const { return scalarBasisFunctionSet().entity(); }

    private:
      template< class Evaluate, class Point, class DofVector >
      void axpy ( const Point &x, const typename Evaluate::Vector &factor, DofVector &dofs ) const
      {
        const std::size_t size = scalarBasisFunctionSet.size();
        std::vector< typename Evaluate::Scalar > scalars( size );
        Evaluate::apply( scalarBasisFunctionSet(), x, scalars );
        for( std::size_t i = 0; i < size; ++i )
        {
          for( int r = 0; r < dimRange; ++r )
          {
            const std::size_t index = r*size + i;
            dofs[ index ] += factor[ r ]*scalars[ i ];
          }
        }
      }

      template< class Evaluate, class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, typename Evaluate::Vector &vectorial )
      {
        const std::size_t size = scalarBasisFunctionSet().size();
        typename Evaluate::Scalar scalar;
        for( int i = 0; i < dimRange; ++i )
        {
          Evaluate::apply( scalarBasisFunctionSet(), x, &(dofs[ i*size ]), scalar );
          vectorial[ i ] = scalar[ 0 ];
        }
      }

      template< class Evaluate, class Point, class VectorArray >
      void evaluateAll ( const Point &x, VectorArray &vectorials ) const
      {
        const std::size_t size = scalarBasisFunctionSet().size();
        std::vector< typename Evaluate::Scalar > scalars( size );
        Evaluate::apply( scalarBasisFunctionSet(), x, scalars );

        typedef typename Evaluate::Vector Vector;

        for( std::size_t i = 0; i < size; ++i )
        {
          for( int r = 0; r < dimRange; ++r )
          {
            const std::size_t index = r*size + i;
            Vector &vectorial = vectorials[ index ];
            vectorial = Vector( typename Vector::value_type( 0 ) );
            vectorial[ r ] = scalars[ i ][ 0 ];
          }
        }
      }

      const ScalarBasisFunctionSetType &scalarBasisFunctionSet () const
      {
        return scalarBasisFunctionSet_;
      }

      ScalarBasisFunctionSetType scalarBasisFunctionSet_;
    };



    // VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::Evaluate
    // --------------------------------------------------------------------

    template< class ScalarBasisFunctionSet, class Range >
    struct VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::Evaluate
    {
      typedef typename ScalarFunctionSpaceType::RangeType Scalar;
      typedef RangeType Vector;

      template< class Point, class DofVector >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, const RangeFieldType *dofs, Scalar &scalar )
      {
        scalarBasisFunctionSet.evaluateAll( x, dofs, scalar );
      }

      template< class Point, class ScalarArray >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, ScalarArray &scalars )
      {
        scalarBasisFunctionSet.evaluateAll( x, scalars );
      }
    };



    // VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::Jacobian
    // --------------------------------------------------------------------

    template< class ScalarBasisFunctionSet, class Range >
    struct VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::Jacobian
    {
      typedef typename ScalarFunctionSpaceType::JacobianRangeType Scalar;
      typedef JacobianRangeType Vector;

      template< class Point, class DofVector >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, const RangeFieldType *dofs, Scalar &scalar )
      {
        scalarBasisFunctionSet.jacobianAll( x, dofs, scalar );
      }

      template< class Point, class ScalarArray >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, ScalarArray &scalars )
      {
        scalarBasisFunctionSet.jacobianAll( x, scalars );
      }
    };



    // VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::Hessian
    // -------------------------------------------------------------------

    template< class ScalarBasisFunctionSet, class Range >
    struct VectorialBasisFunctionSet< ScalarBasisFunctionSet, Range >::Hessian
    {
      typedef typename ScalarFunctionSpaceType::HessianRangeType Scalar;
      typedef HessianRangeType Vector;

      template< class Point, class DofVector >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, const RangeFieldType *dofs, Scalar &scalar )
      {
        scalarBasisFunctionSet.hessianAll( x, dofs, scalar );
      }

      template< class Point, class ScalarArray >
      static void apply ( const ScalarBasisFunctionSetType &scalarBasisFunctionSet,
                          const Point &x, ScalarArray &scalars )
      {
        scalarBasisFunctionSet.hessianAll( x, scalars );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BASISFUNCTIONSET_VECTORIAL_HH
