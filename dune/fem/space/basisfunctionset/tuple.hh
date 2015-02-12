#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TUPLE_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_TUPLE_HH

#include <algorithm>
#include <array>
#include <tuple>
#include <vector>

#include <dune/common/forloop.hh>
#include <dune/common/std/utility.hh>
#include <dune/common/tuples.hh>

#include <dune/geometry/type.hh>

#include <dune/fem/common/utility.hh>

#include <dune/fem/space/basisfunctionset/basisfunctionset.hh>
#include <dune/fem/space/combinedspace/subobjects.hh>


namespace Dune
{

  namespace Fem
  {

    // TupleBasisFunctionSet
    // ---------------------

    template< class ... BasisFunctionSets >
    class TupleBasisFunctionSet
    {
      typedef TupleBasisFunctionSet< BasisFunctionSets ... > ThisType;

      // Functor Helper structs
      template< int > struct EvaluateAll;
      template< int > struct JacobianAll;
      template< int > struct HessianAll;
      template< int > struct EvaluateAllRanges;
      template< int > struct JacobianAllRanges;
      template< int > struct HessianAllRanges;
      template< int > struct Axpy;

      // helper class to compute static overall size and offsets
      template< int ... i >
      struct Indices
      {
        template< int j >
        static constexpr int size () { return Dune::tuple_element< j, Dune::tuple< std::integral_constant< int, i > ... > >::type::value; }

        template< int ... j >
        static constexpr Std::integer_sequence< int, size< j >() ... > sizes ( Std::integer_sequence< int, j ... > )
        {
          return Std::integer_sequence< int, size< j >() ... >();
        }

        template< int j >
        static constexpr int offset ()
        {
          return mysum( sizes( Std::make_integer_sequence< int, j >() ) );
        }

      private:
        template< int ... j >
        static constexpr int mysum ( Std::integer_sequence< int, j ... > )
        {
          return Std::sum( j ... );
        }

        static constexpr int mysum ( Std::integer_sequence< int > )
        {
          return 0;
        }
      };

      //! helper class to compute static rangeoffsets
      typedef Indices< BasisFunctionSets::FunctionSpaceType::dimRange ... > RangeIndices;

      //! type of tuple
      typedef Dune::tuple< BasisFunctionSets ... > BasisFunctionSetTupleType;

      static_assert( Std::are_all_same< typename BasisFunctionSets::DomainType ... >::value, 
          "TupleBasisFunctionSet needs common DomainType" );
      //! get type of first function space, to obtain dimDomain and Field Types
      typedef typename Dune::tuple_element< 0, BasisFunctionSetTupleType >::type::FunctionSpaceType ContainedFunctionSpaceType;

      //! number of BasisFunctionsSets in the tuple
      static const int setSize = sizeof ... ( BasisFunctionSets ) -1;

    public:
      //! size of domian space
      static const int dimDomain = ContainedFunctionSpaceType::dimDomain;

      //! size of range space
      static constexpr const int dimRange = Std::sum( static_cast< int >( BasisFunctionSets::FunctionSpaceType::dimRange ) ... );

      //! type of analytical combiend function space
      typedef FunctionSpace< typename ContainedFunctionSpaceType::DomainFieldType,
                             typename ContainedFunctionSpaceType::RangeFieldType, dimDomain, dimRange > FunctionSpaceType;

      //! type of Domain Vector
      typedef typename FunctionSpaceType::DomainType DomainType;

      //! type of Range Vector
      typedef typename FunctionSpaceType::RangeType RangeType;

      //! type of Range Vector field
      typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

      //! type of Jacobian Vector/Matrix
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

      //! type of Hessian Matrix
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      static_assert( Std::are_all_same< typename BasisFunctionSets::EntityType ... >::value, 
          "TupleBasisFunctionSet needs common EntityType" );
      //! type of Entity the basis function set is initialized on
      typedef typename Dune::tuple_element< 0, BasisFunctionSetTupleType >::type::EntityType EntityType;

      static_assert( Std::are_all_same< typename BasisFunctionSets::ReferenceElementType ... >::value, 
          "TupleBasisFunctionSet needs ReferenceElementType" );
      //! type of reference element for this BasisFunctionSet
      typedef typename Dune::tuple_element< 0, BasisFunctionSetTupleType >::type::ReferenceElementType ReferenceElementType;


      // default constructor, needed by localmatrix
      TupleBasisFunctionSet () {}

      // constructor taking a pack of basisFunctionSets
      TupleBasisFunctionSet ( const BasisFunctionSets & ... basisFunctionSets )
        : basisFunctionSetTuple_( basisFunctionSets ... )
      {}

      // constructor taking a tuple of basisfunction sets
      TupleBasisFunctionSet ( const BasisFunctionSetTupleType &basisFunctionSetTuple )
        : basisFunctionSetTuple_( basisFunctionSetTuple )
      {}

      // Basis Function Set Interface Methods
      // ------------------------------------

      //! \brief return order of basis function set, maximal order in the tupleset
      int order () const
      {
        return order( Std::index_sequence_for< BasisFunctionSets ... >() );
      }

      //! \brief return size of basis function set
      std::size_t size () const
      {
        return size( Std::index_sequence_for< BasisFunctionSets ... >() );
      }

      //! \copydoc BasisFunctionSet::type
      Dune::GeometryType type () const
      {
        return std::get< 0 >( basisFunctionSetTuple_ ).type();
      }

      //! \copydoc BasisFunctionSet::entity
      const EntityType &entity () const
      {
        return std::get< 0 >( basisFunctionSetTuple_ ).entity();
      }

      //! \copydoc BasisFunctionSet::entity
      const ReferenceElementType &referenceElement () const
      {
        return std::get< 0 >( basisFunctionSetTuple_ ).referenceElement();
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, dofs, value )
      template< class Point, class DofVector >
      void evaluateAll ( const Point &x, const DofVector &dofs, RangeType &value ) const
      {
        std::size_t offset = 0;
        ForLoop< EvaluateAll, 0, setSize >::apply( offset, x, dofs, value, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, values )
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        assert( values.size() >= size() );
        std::size_t offset = 0;
        ForLoop< EvaluateAllRanges, 0, setSize >::apply( offset, x, values, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::evaluateAll( quad, dofs, ranges )
      template< class QuadratureType, class DofVector, class RangeArray >
      void evaluateAll ( const QuadratureType &quad, const DofVector &dofs, RangeArray &ranges ) const
      {
        const int nop = quad.nop();
        for( int qp = 0; qp < nop; ++qp )
          evaluateAll( quad[ qp ], dofs, ranges[ qp ] );
      }

      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobian )
      template< class Point, class DofVector >
      void jacobianAll ( const Point &x, const DofVector &dofs, JacobianRangeType &jacobian ) const
      {
        std::size_t offset = 0;
        ForLoop< JacobianAll, 0, setSize >::apply( offset, x, dofs, jacobian, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobians )
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        assert( jacobians.size() >= size() );
        std::size_t offset = 0;
        ForLoop< JacobianAllRanges, 0, setSize >::apply( offset, x, jacobians, basisFunctionSetTuple_ );
      }

      //! \brief evaluate the jacobian of all basis functions and store the result in the jacobians array
      template< class QuadratureType, class DofVector, class JacobianArray >
      void jacobianAll ( const QuadratureType &quad, const DofVector &dofs, JacobianArray &jacobians ) const
      {
        const int nop = quad.nop();
        for( int qp = 0; qp < nop; ++qp )
          jacobianAll( quad[ qp ], dofs, jacobians[ qp ] );
      }

      //! \copydoc BasisFunctionSet::hessianAll( x, dofs, hessian )
      template< class Point, class DofVector >
      void hessianAll ( const Point &x, const DofVector &dofs, HessianRangeType &hessian ) const
      {
        std::size_t offset = 0;
        ForLoop< HessianAll, 0, setSize >::apply( offset, x, dofs, hessian, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::hessianAll( x, hessians )
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        assert( hessians.size() >= size() );
        std::size_t offset = 0;
        ForLoop< HessianAllRanges, 0, setSize >::apply( offset, x, hessians, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::axpy( quad, values, dofs )
      template< class QuadratureType, class Vector, class DofVector >
      void axpy ( const QuadratureType &quad, const Vector &values, DofVector &dofs ) const
      {
        // call axpy method for each entry of the given vector, e.g. rangeVector or jacobianVector
        const unsigned int nop = quad.nop();
        for( unsigned int qp = 0; qp < nop; ++qp )
          axpy( quad[ qp ], values[ qp ], dofs );
      }

      //! \copydoc BasisFunctionSet::axpy( quad, valuesA, valuesB, dofs )
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

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, DofVector &dofs ) const
      {
        std::size_t offset = 0;
        ForLoop< Axpy, 0, setSize >::apply( offset, x, valueFactor, dofs, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::axpy( x, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        std::size_t offset = 0;
        ForLoop< Axpy, 0, setSize >::apply( offset, x, jacobianFactor, dofs, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        std::size_t offset = 0;
        ForLoop< Axpy, 0, setSize >::apply( offset, x, valueFactor, jacobianFactor, dofs, basisFunctionSetTuple_ );
      }

    protected:
      // unroll index sequence and take maxmial order
      template< std::size_t ... i >
      int order ( Std::index_sequence< i ... > ) const
      {
        return Std::max( std::get< i >( basisFunctionSetTuple_ ).order() ... );
      }

      // unroll index sequence and sum up sizes
      template< std::size_t ... i >
      std::size_t size ( Std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( basisFunctionSetTuple_ ).size() ... );
      }

    private:
      BasisFunctionSetTupleType basisFunctionSetTuple_;
    };



    // EvaluateAll
    // -----------

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    EvaluateAll
    {
      static const int rangeOffset = RangeIndices::template offset< i >();

      template< class Point, class DofVector, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, const DofVector &dofVector, RangeType &values, const Tuple &tuple )
      {
        // get size of this basisFunctionSet
        std::size_t size = std::get< i >( tuple ).size();
        typedef typename DofVector::value_type DofType;
        typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::RangeType thisRange;

        SubDofVector< const DofVector, const DofType > subDofVector( dofVector, size, offset );
        // evaluateAll for this BasisFunctionSet, with View on DofVector
        std::get< i >( tuple ).evaluateAll( x, subDofVector, thisRange );
        std::copy( thisRange.begin(), thisRange.end(), values.begin() + rangeOffset );
        offset += size;
      }
    };


    // JacobianAll
    // -----------

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    JacobianAll
    {
      static const int rangeOffset = RangeIndices::template offset< i >();

      template< class Point, class DofVector, class Tuple >
      static void apply ( std::size_t &offset, const Point &x,
                          const DofVector &dofVector, JacobianRangeType &values, const Tuple &tuple )
      {
        // get size of this basisFunctionSet
        std::size_t size = std::get< i >( tuple ).size();
        typedef typename DofVector::value_type DofType;
        typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::JacobianRangeType thisJacobian;

        SubDofVector< const DofVector, const DofType > subDofVector( dofVector, size, offset );
        // jacobianAll for this BasisFunctionSet, with View on DofVector
        std::get< i >( tuple ).jacobianAll( x, subDofVector, thisJacobian );
        std::copy( thisJacobian.begin(), thisJacobian.end(), values.begin() + rangeOffset );
        offset += size;
      }
    };


    // HessianAll
    // ----------

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    HessianAll
    {
      static const int rangeOffset = RangeIndices::template offset< i >();

      template< class Point, class DofVector, class Tuple >
      static void apply ( std::size_t &offset, const Point &x,
                          const DofVector &dofVector, HessianRangeType &values, const Tuple &tuple )
      {
        // get size of this basisFunctionSet
        std::size_t size = std::get< i >( tuple ).size();
        typedef typename DofVector::value_type DofType;
        typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::HessianRangeType thisHessian;

        SubDofVector< const DofVector, const DofType > subDofVector( dofVector, size, offset );
        // hessianAll for this BasisFunctionSet, with View on DofVector
        std::get< i >( tuple ).hessianAll( x, subDofVector, thisHessian );
        std::copy( thisHessian.begin(), thisHessian.end(), values.begin() + rangeOffset );
        offset += size;
      }
    };


    // EvaluateAllRanges
    // -----------------

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    EvaluateAllRanges
    {
      typedef typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::RangeType ThisRangeType;
      static const int thisDimRange = Dune::tuple_element< i, BasisFunctionSetTupleType >::type::FunctionSpaceType::dimRange;
      static const int rangeOffset = RangeIndices::template offset< i >();

      template< class Point, class RangeArray, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, RangeArray &values, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        std::vector< ThisRangeType > thisValues( size );

        // misuse SubDofVector
        std::get< i >( tuple ).evaluateAll( x, thisValues );
        for( std::size_t j = 0; j < size; ++j )
        {
          values[ j + offset ] = RangeType( 0.0 );
          for( int r = 0; r < thisDimRange; ++r )
            values[ j + offset ][ r + rangeOffset ] = thisValues[ j ][ r ];
        }

        offset += size;
      }
    };


    // JacobianAllRanges
    // -----------------

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    JacobianAllRanges
    {
      typedef typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::JacobianRangeType ThisJacobianRangeType;
      static const int thisDimRange = Dune::tuple_element< i, BasisFunctionSetTupleType >::type::FunctionSpaceType::dimRange;
      static const int rangeOffset = RangeIndices::template offset< i >();

      template< class Point, class RangeArray, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, RangeArray &values, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();

        std::vector< ThisJacobianRangeType > thisValues( size );
        std::get< i >( tuple ).jacobianAll( x, thisValues );

        for( std::size_t j = 0; j < size; ++j )
        {
          values[ j + offset ] = JacobianRangeType( RangeFieldType( 0.0 ) );
          for( int r = 0; r < thisDimRange; ++r )
            values[ j + offset ][ r + rangeOffset ] = thisValues[ j ][ r ];
        }

        offset += size;
      }
    };


    // HessianAllRanges
    // ----------------

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    HessianAllRanges
    {
      typedef typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::HessianRangeType ThisHessianRangeType;
      static const int thisDimRange = Dune::tuple_element< i, BasisFunctionSetTupleType >::type::FunctionSpaceType::dimRange;
      static const int rangeOffset = RangeIndices::template offset< i >();

      template< class Point, class RangeArray, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, RangeArray &values, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();

        std::vector< ThisHessianRangeType > thisValues( size );
        std::get< i >( tuple ).hessianAll( x, thisValues );

        for( std::size_t j = 0; j < size; ++j )
        {
          values[ j + offset ] = HessianRangeType( RangeFieldType( 0.0 ) );
          for( int r = 0; r < thisDimRange; ++r )
            values[ j + offset ][ r + rangeOffset ] = thisValues[ j ][ r ];
        }

        offset += size;
      }
    };


    // Axpy
    // ----

    template< class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< BasisFunctionSets ... >::
    Axpy
    {
      typedef typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::RangeType ThisRangeType;
      typedef typename Dune::tuple_element< i, BasisFunctionSetTupleType >::type::JacobianRangeType ThisJacobianRangeType;

      static const int rangeOffset = RangeIndices::template offset< i >();

      typedef SubObject< const RangeType, const ThisRangeType, rangeOffset > SubRangeType;
      typedef SubObject< const JacobianRangeType, const ThisJacobianRangeType, rangeOffset > SubJacobianRangeType;

      // axpy with range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, const RangeType &factor, DofVector &dofVector, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubRangeType subFactor( factor );
        SubDofVector< DofVector, typename DofVector::value_type > subDofVector( dofVector, size, offset );
        std::get< i >( tuple ).axpy( x, (ThisRangeType) subFactor, subDofVector );
        offset += size;
      }

      // axpy with jacobian range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, const JacobianRangeType &factor, DofVector &dofVector, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubJacobianRangeType subFactor( factor );
        SubDofVector< DofVector, typename DofVector::value_type > subDofVector( dofVector, size, offset );
        std::get< i >( tuple ).axpy( x, (ThisJacobianRangeType) subFactor, subDofVector );
        offset += size;
      }

      // axpy with range and jacobian range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( std::size_t &offset, const Point &x, const RangeType &rangeFactor, const JacobianRangeType &jacobianFactor, DofVector &dofVector, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubRangeType subRangeFactor( rangeFactor );
        SubJacobianRangeType subJacobianFactor( jacobianFactor );
        SubDofVector< DofVector, typename DofVector::value_type > subDofVector( dofVector, size, offset );
        std::get< i >( tuple ).axpy( x, (ThisRangeType) subRangeFactor, (ThisJacobianRangeType) subJacobianFactor, subDofVector );
        offset += size;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TUPLE_HH
