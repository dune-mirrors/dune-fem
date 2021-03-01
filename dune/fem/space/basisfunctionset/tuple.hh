#ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TUPLE_HH
#define DUNE_FEM_SPACE_BASISFUNCTIONSET_TUPLE_HH

#include <algorithm>
#include <array>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/fem/common/forloop.hh>

#include <dune/geometry/type.hh>

#include <dune/fem/common/utility.hh>

#include <dune/fem/space/basisfunctionset/basisfunctionset.hh>
#include <dune/fem/space/combinedspace/subobjects.hh>
#include <dune/fem/storage/subvector.hh>

namespace Dune
{

  namespace Fem
  {

    // TupleBasisFunctionSet
    // ---------------------
    //
    // TupleDiscreteFunctionSpace combination operation
    // describes the way in which the spaces have been combined, i.e.
    // Product:    V = V_1 x V_2 x ...
    // Summation:  V = V_1 + V_2 + ...

    struct TupleSpaceProduct {};
    struct TupleSpaceSummation {};

    template< class CombineOp, class ... BasisFunctionSets >
    class TupleBasisFunctionSet
    {
      typedef TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... > ThisType;

      // Functor Helper structs
      template< int > struct ComputeOffset;
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
        static constexpr int size () { return std::tuple_element< j, std::tuple< std::integral_constant< int, i > ... > >::type::value; }

        template< int ... j >
        static constexpr std::integer_sequence< int, size< j >() ... > sizes ( std::integer_sequence< int, j ... > )
        {
          return std::integer_sequence< int, size< j >() ... >();
        }

        template< int j >
        static constexpr int offset ()
        {
          return mysum( sizes( std::make_integer_sequence< int, j >() ) );
        }

      private:
        template< int ... j >
        static constexpr int mysum ( std::integer_sequence< int, j ... > )
        {
          return Std::sum( j ... );
        }

        static constexpr int mysum ( std::integer_sequence< int > )
        {
          return 0;
        }
      };

      //! type of tuple
      typedef std::tuple< BasisFunctionSets ... > BasisFunctionSetTupleType;

      static_assert( Std::are_all_same< typename BasisFunctionSets::DomainType ... >::value,
          "TupleBasisFunctionSet needs common DomainType" );
      //! get type of first function space, to obtain dimDomain and Field Types
      typedef typename std::tuple_element< 0, BasisFunctionSetTupleType >::type::FunctionSpaceType ContainedFunctionSpaceType;

      //! size of domain space
      static const int dimDomain = ContainedFunctionSpaceType::dimDomain;

      //! number of BasisFunctionsSets in the tuple
      static const int setSize = sizeof ... ( BasisFunctionSets );
      static const int setIterationSize = sizeof ... ( BasisFunctionSets )-1;

      //! type of offset array
      typedef std::array< std::size_t, setSize + 1 > OffsetType;

      // storage for i-th sub basisfunction
      struct Storage
      {
        typedef std::tuple< std::vector< typename BasisFunctionSets::RangeType > ... > RangeStorageTupleType;
        typedef std::tuple< std::vector< typename BasisFunctionSets::JacobianRangeType > ... > JacobianStorageTupleType;
        typedef std::tuple< std::vector< typename BasisFunctionSets::HessianRangeType > ... > HessianStorageTupleType;

        Storage () {}

        template< int i >
        typename std::tuple_element< i, RangeStorageTupleType >::type
        & rangeStorage() const { return std::get< i >( rangeStorage_ ); }

        template< int i >
        typename std::tuple_element< i, JacobianStorageTupleType >::type
        & jacobianStorage() const { return std::get< i >( jacobianStorage_ ); }

        template< int i >
        typename std::tuple_element< i, HessianStorageTupleType >::type
        & hessianStorage() const { return std::get< i >( hessianStorage_ ); }

      private:
        mutable RangeStorageTupleType rangeStorage_;
        mutable JacobianStorageTupleType jacobianStorage_;
        mutable HessianStorageTupleType hessianStorage_;
      };


    public:
      //! helper class to compute static rangeoffsets
      typedef Indices< BasisFunctionSets::FunctionSpaceType::dimRange ... > RangeIndices;

    protected:
      template <int dummy, class CombOp>
      struct CombinationTraits;

      template <int dummy>
      struct CombinationTraits< dummy, TupleSpaceProduct >
      {
        //! size of range space for product space
        static constexpr const int dimRange = Std::sum( static_cast< int >( BasisFunctionSets::FunctionSpaceType::dimRange ) ... );

        //! type of analytical combined product function space
        typedef  FunctionSpace< typename ContainedFunctionSpaceType::DomainFieldType,
                                typename ContainedFunctionSpaceType::RangeFieldType,
                                dimDomain, dimRange > FunctionSpaceType;

        template <class ThisRange, class Range>
        static void apply(const int rangeOffset, const ThisRange& thisRange, Range& values)
        {
          // scatter result to correct place in larger result vector
          std::copy( thisRange.begin(), thisRange.end(), values.begin() + rangeOffset );
        }

        template< int i >
        static constexpr int rangeOffset ()
        {
          return RangeIndices::template offset< i >();
        }
      };

      template <int dummy>
      struct CombinationTraits< dummy, TupleSpaceSummation >
      {
        //! type of analytical is same as contained space since its a sum
        typedef ContainedFunctionSpaceType  FunctionSpaceType;

        template <class ThisRange, class Range>
        static void apply(const int rangeOffset, const ThisRange& thisRange, Range& values)
        {
          // sum results
          values += thisRange;
        }

        template< int j >
        static constexpr int rangeOffset ()
        {
          return 0;
        }
      };

      // type of accumulation of result after loop over basis sets
      typedef CombinationTraits< 0, CombineOp > CombinationType;

    public:
      // export type of i-th subbasisfunction set
      template< int i >
      struct SubBasisFunctionSet
      {
        typedef typename std::tuple_element< i, BasisFunctionSetTupleType >::type type;
      };

      //! type of analytical combined function space
      typedef typename CombinationType :: FunctionSpaceType  FunctionSpaceType;

      //! size of domain space
      static const int dimRange = FunctionSpaceType::dimRange;

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
      typedef typename std::tuple_element< 0, BasisFunctionSetTupleType >::type::EntityType EntityType;

      static_assert( Std::are_all_same< typename BasisFunctionSets::ReferenceElementType ... >::value,
          "TupleBasisFunctionSet needs ReferenceElementType" );
      //! type of reference element for this BasisFunctionSet
      typedef typename std::tuple_element< 0, BasisFunctionSetTupleType >::type::ReferenceElementType ReferenceElementType;

      // default constructor, needed by localmatrix
      TupleBasisFunctionSet () : offset_( {{ 0 }} ) {}

      // constructor taking a pack of basisFunctionSets
      TupleBasisFunctionSet ( const BasisFunctionSets & ... basisFunctionSets )
        : basisFunctionSetTuple_( basisFunctionSets ... ),
          offset_()
      {
        offset_[ 0 ] = 0;
        Fem::ForLoop< ComputeOffset, 0, setIterationSize >::apply( offset_, basisFunctionSetTuple_ );
      }

      // constructor taking a tuple of basisfunction sets
      TupleBasisFunctionSet ( const BasisFunctionSetTupleType &basisFunctionSetTuple )
        : basisFunctionSetTuple_( basisFunctionSetTuple ),
          offset_()
      {
        offset_[ 0 ] = 0;
        Fem::ForLoop< ComputeOffset, 0, setIterationSize >::apply( offset_, basisFunctionSetTuple_ );
      }

      // Basis Function Set Interface Methods
      // ------------------------------------

      //! \brief return order of basis function set, maximal order in the tupleset
      int order () const
      {
        return order( std::index_sequence_for< BasisFunctionSets ... >() );
      }

      //! \brief return size of basis function set
      std::size_t size () const
      {
        return size( std::index_sequence_for< BasisFunctionSets ... >() );
      }

      //! \copydoc BasisFunctionSet::type
      Dune::GeometryType type () const
      {
        return std::get< 0 >( basisFunctionSetTuple_ ).type();
      }

      //! \copydoc BasisFunctionSet::valid
      bool valid () const
      {
        return std::get< 0 >( basisFunctionSetTuple_ ).valid();
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
        Fem::ForLoop< EvaluateAll, 0, setIterationSize >::apply( x, dofs, value, offset_, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::evaluateAll( x, values )
      template< class Point, class RangeArray >
      void evaluateAll ( const Point &x, RangeArray &values ) const
      {
        assert( values.size() >= size() );
        Fem::ForLoop< EvaluateAllRanges, 0, setIterationSize >::apply( x, values, storage_, offset_, basisFunctionSetTuple_ );
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
        Fem::ForLoop< JacobianAll, 0, setIterationSize >::apply( x, dofs, jacobian, offset_, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::jacobianAll( x, dofs, jacobians )
      template< class Point, class JacobianRangeArray >
      void jacobianAll ( const Point &x, JacobianRangeArray &jacobians ) const
      {
        assert( jacobians.size() >= size() );
        Fem::ForLoop< JacobianAllRanges, 0, setIterationSize >::apply( x, jacobians, storage_, offset_, basisFunctionSetTuple_ );
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
        Fem::ForLoop< HessianAll, 0, setIterationSize >::apply( x, dofs, hessian, offset_, basisFunctionSetTuple_ );
      }

      //! \brief evaluate the hessian of all basis functions and store the result in the hessians array
      template< class QuadratureType, class DofVector, class HessianArray >
      void hessianAll ( const QuadratureType &quad, const DofVector &dofs, HessianArray &hessians ) const
      {
        const int nop = quad.nop();
        for( int qp = 0; qp < nop; ++qp )
          hessianAll( quad[ qp ], dofs, hessians[ qp ] );
      }

      //! \copydoc BasisFunctionSet::hessianAll( x, hessians )
      template< class Point, class HessianRangeArray >
      void hessianAll ( const Point &x, HessianRangeArray &hessians ) const
      {
        assert( hessians.size() >= size() );
        Fem::ForLoop< HessianAllRanges, 0, setIterationSize >::apply( x, hessians, storage_, offset_, basisFunctionSetTuple_ );
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
        Fem::ForLoop< Axpy, 0, setIterationSize >::apply( x, valueFactor, dofs, offset_, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::axpy( x, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        Fem::ForLoop< Axpy, 0, setIterationSize >::apply( x, jacobianFactor, dofs, offset_, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::axpy( x, hessianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const HessianRangeType &hessianFactor, DofVector &dofs ) const
      {
        Fem::ForLoop< Axpy, 0, setIterationSize >::apply( x, hessianFactor, dofs, offset_, basisFunctionSetTuple_ );
      }

      //! \copydoc BasisFunctionSet::axpy( x, valueFactor, jacobianFactor, dofs )
      template< class Point, class DofVector >
      void axpy ( const Point &x, const RangeType &valueFactor, const JacobianRangeType &jacobianFactor, DofVector &dofs ) const
      {
        Fem::ForLoop< Axpy, 0, setIterationSize >::apply( x, valueFactor, jacobianFactor, dofs, offset_, basisFunctionSetTuple_ );
      }

      /***** NON Interface methods ****/

      //! return i-th subbasisfunctionSet
      template< int i >
      const typename SubBasisFunctionSet< i >::type& subBasisFunctionSet() const
      {
        return std::get< i >( basisFunctionSetTuple_ );
      }

      //! return offset of the i-th subbasisfunctionSet in the whole set
      std::size_t offset ( int i ) const
      {
        return offset_[ i ];
      }

      //! return number of subBasisFunctionSets
      static int numSubBasisFunctionSets ()
      {
        return setSize;
      }

    protected:
      // unroll index sequence and take maximal order
      template< std::size_t ... i >
      int order ( std::index_sequence< i ... > ) const
      {
        return Std::max( std::get< i >( basisFunctionSetTuple_ ).order() ... );
      }

      // unroll index sequence and sum up sizes
      template< std::size_t ... i >
      std::size_t size ( std::index_sequence< i ... > ) const
      {
        return Std::sum( std::get< i >( basisFunctionSetTuple_ ).size() ... );
      }

    private:
      BasisFunctionSetTupleType basisFunctionSetTuple_;
      OffsetType offset_;

      Storage storage_;
    };



    // ComputeOffset
    // -------------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    ComputeOffset
    {
      template< class Tuple >
      static void apply ( OffsetType &offset, const Tuple &tuple )
      {
        offset[ i + 1 ] = offset[ i ] + std::get< i >( tuple ).size();
      }
    };


    // EvaluateAll
    // -----------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    EvaluateAll
    {
      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const DofVector &dofVector, RangeType &values, const OffsetType &offset, const Tuple &tuple )
      {
        // evaluateAll for this BasisFunctionSet, with View on DofVector
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< const DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        typename std::tuple_element< i, BasisFunctionSetTupleType >::type::RangeType thisRange;
        std::get< i >( tuple ).evaluateAll( x, subDofVector, thisRange );

        // distribute result to values
        CombinationType::apply( rangeOffset, thisRange, values );
      }
    };


    // JacobianAll
    // -----------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    JacobianAll
    {
      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const DofVector &dofVector, JacobianRangeType &values, const OffsetType &offset, const Tuple &tuple )
      {
        // jacobianAll for this BasisFunctionSet, with View on DofVector
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< const DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        typename std::tuple_element< i, BasisFunctionSetTupleType >::type::JacobianRangeType thisJacobian;
        std::get< i >( tuple ).jacobianAll( x, subDofVector, thisJacobian );

        // distribute result to values
        CombinationType::apply( rangeOffset, thisJacobian, values );
      }
    };


    // HessianAll
    // ----------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    HessianAll
    {
      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const DofVector &dofVector, HessianRangeType &values, const OffsetType &offset, const Tuple &tuple )
      {
        // hessianAll for this BasisFunctionSet, with View on DofVector
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< const DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        typename std::tuple_element< i, BasisFunctionSetTupleType >::type::HessianRangeType thisHessian;
        std::get< i >( tuple ).hessianAll( x, subDofVector, thisHessian );

        // distribute result to values
        CombinationType::apply( rangeOffset, thisHessian, values );
      }
    };


    // EvaluateAllRanges
    // -----------------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    EvaluateAllRanges
    {
      typedef typename std::tuple_element< i, typename Storage::RangeStorageTupleType >::type ThisStorage;
      static const int thisDimRange = std::tuple_element< i, BasisFunctionSetTupleType >::type::FunctionSpaceType::dimRange;

      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      template< class Point, class RangeArray, class Tuple >
      static void apply ( const Point &x, RangeArray &values, const Storage &s, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        ThisStorage &thisStorage = s.template rangeStorage< i >();
        thisStorage.resize( size );

        std::get< i >( tuple ).evaluateAll( x, thisStorage );

        for( std::size_t j = 0; j < size; ++j )
        {
          values[ j + offset[ i ] ] = RangeType( 0.0 );
          for( int r = 0; r < thisDimRange; ++r )
            values[ j + offset[ i ] ][ r + rangeOffset ] = thisStorage[ j ][ r ];
        }
      }
    };


    // JacobianAllRanges
    // -----------------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    JacobianAllRanges
    {
      typedef typename std::tuple_element< i, typename Storage::JacobianStorageTupleType >::type ThisStorage;
      static const int thisDimRange = std::tuple_element< i, BasisFunctionSetTupleType >::type::FunctionSpaceType::dimRange;

      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      template< class Point, class RangeArray, class Tuple >
      static void apply ( const Point &x, RangeArray &values, const Storage &s, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        ThisStorage &thisStorage = s.template jacobianStorage< i >();
        thisStorage.resize( size );

        std::get< i >( tuple ).jacobianAll( x, thisStorage );

        for( std::size_t j = 0; j < size; ++j )
        {
          values[ j + offset[ i ] ] = JacobianRangeType( RangeFieldType( 0.0 ) );
          for( int r = 0; r < thisDimRange; ++r )
            values[ j + offset[ i ] ][ r + rangeOffset ] = thisStorage[ j ][ r ];
        }
      }
    };


    // HessianAllRanges
    // ----------------

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    HessianAllRanges
    {
      typedef typename std::tuple_element< i, typename Storage::HessianStorageTupleType >::type ThisStorage;
      static const int thisDimRange = std::tuple_element< i, BasisFunctionSetTupleType >::type::FunctionSpaceType::dimRange;

      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      template< class Point, class RangeArray, class Tuple >
      static void apply ( const Point &x, RangeArray &values, const Storage &s, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        ThisStorage &thisStorage = s.template hessianStorage< i >();
        thisStorage.resize( size );

        std::get< i >( tuple ).hessianAll( x, thisStorage );

        for( std::size_t j = 0; j < size; ++j )
        {
          values[ j + offset[ i ] ] = typename HessianRangeType::value_type( 0.0 );
          for( int r = 0; r < thisDimRange; ++r )
            values[ j + offset[ i ] ][ r + rangeOffset ] = thisStorage[ j ][ r ];
        }
      }
    };


    // Axpy
    // ----

    template< class CombineOp, class ... BasisFunctionSets >
    template< int i >
    struct TupleBasisFunctionSet< CombineOp, BasisFunctionSets ... >::
    Axpy
    {
      typedef typename std::tuple_element< i, BasisFunctionSetTupleType >::type::RangeType ThisRangeType;
      typedef typename std::tuple_element< i, BasisFunctionSetTupleType >::type::JacobianRangeType ThisJacobianRangeType;
      typedef typename std::tuple_element< i, BasisFunctionSetTupleType >::type::HessianRangeType ThisHessianRangeType;

      // only needed for product spaces
      static const int rangeOffset = CombinationType :: template rangeOffset< i >();

      typedef SubObject< const RangeType, const ThisRangeType, rangeOffset > SubRangeType;
      typedef SubObject< const JacobianRangeType, const ThisJacobianRangeType, rangeOffset > SubJacobianRangeType;
      typedef SubObject< const HessianRangeType, const ThisHessianRangeType, rangeOffset > SubHessianRangeType;

      // axpy with range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const RangeType &factor, DofVector &dofVector, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        SubRangeType subFactor( factor );
        std::get< i >( tuple ).axpy( x, (ThisRangeType) subFactor, subDofVector );
      }

      // axpy with jacobian range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const JacobianRangeType &factor, DofVector &dofVector, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        SubJacobianRangeType subFactor( factor );
        std::get< i >( tuple ).axpy( x, (ThisJacobianRangeType) subFactor, subDofVector );
      }

      // axpy with hessian range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const HessianRangeType &factor, DofVector &dofVector, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        SubHessianRangeType subFactor( factor );
        std::get< i >( tuple ).axpy( x, (ThisHessianRangeType) subFactor, subDofVector );
      }

      // axpy with range and jacobian range type
      template< class Point, class DofVector, class Tuple >
      static void apply ( const Point &x, const RangeType &rangeFactor, const JacobianRangeType &jacobianFactor, DofVector &dofVector, const OffsetType &offset, const Tuple &tuple )
      {
        std::size_t size = std::get< i >( tuple ).size();
        SubVector< DofVector, OffsetSubMapper > subDofVector( dofVector, OffsetSubMapper( size, offset[ i ] ) );
        SubRangeType subRangeFactor( rangeFactor );
        SubJacobianRangeType subJacobianFactor( jacobianFactor );
        std::get< i >( tuple ).axpy( x, (ThisRangeType) subRangeFactor, (ThisJacobianRangeType) subJacobianFactor, subDofVector );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BASISFUNCTIONSET_TUPLE_HH
