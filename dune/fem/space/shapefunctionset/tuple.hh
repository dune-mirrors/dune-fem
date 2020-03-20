#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_TUPLE_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_TUPLE_HH

#include <tuple>

#include <dune/geometry/type.hh>

#include <dune/fem/common/forloop.hh>
#include <dune/common/tupleutility.hh>
#include <dune/fem/common/utility.hh>

#include <dune/fem/space/shapefunctionset/vectorial.hh>


namespace Dune
{

  namespace Fem
  {

    // TupleShapeFunctionSet
    // ---------------------

    template< class ... ShapeFunctionSets >
    class TupleShapeFunctionSet
    {
      typedef TupleShapeFunctionSet< ShapeFunctionSets ... > ThisType;

      template< int ... I >
      struct RangeOffsets
      {
        typedef std::tuple< std::integral_constant< int, I > ... > RangeSizeTuple;

        template< int j >
        static constexpr int size () { return std::tuple_element< j, RangeSizeTuple >::type::value; }

        template< int ... j >
        static constexpr std::integer_sequence< int, size< j >() ... > sizes ( std::integer_sequence< int, j ... > )
        {
          return std::integer_sequence< int, size< j >() ... >();
        }

        template< int i >
        static constexpr int offset ()
        {
          return sum( sizes( std::make_integer_sequence< int, i >() ) );
        }

      private:
        template< int ... j >
        static constexpr int sum ( std::integer_sequence< int, j ... > )
        {
          return Std::sum( j ... );
        }

        static constexpr int sum ( std::integer_sequence< int > ) { return 0; }
      };

      typedef std::array< std::size_t, sizeof ... ( ShapeFunctionSets ) +1 > Offset;

      template< int I > struct Offsets;
      template< class Functor, class Value, int I > struct FunctorWrapper;

      template< int I > struct EvaluateEach;
      template< int I > struct JacobianEach;
      template< int I > struct HessianEach;

      static const std::size_t dimRange = Std::sum( static_cast< int >( ShapeFunctionSets::FunctionSpaceType::dimRange ) ... );

    public:
      template< std::size_t i >
      using SubShapeFunctionSetType = std::tuple_element_t< i, std::tuple< ShapeFunctionSets... > >;

      typedef typename ToNewDimRangeFunctionSpace< typename SubShapeFunctionSetType< 0 >::FunctionSpaceType, dimRange >::Type FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      TupleShapeFunctionSet () : offset_( {{ 0 }} ) {
      }

      TupleShapeFunctionSet ( GeometryType type )
        : shapeFunctionSetTuple_( makeGeometryTypeTuple( type, std::index_sequence_for< ShapeFunctionSets ... >() ) )
      {
        offset_[ 0 ] = 0;
        Fem::ForLoop< Offsets, 0, sizeof ... ( ShapeFunctionSets ) -1 >::apply( shapeFunctionSetTuple_, offset_ );
      }

      template< class ... Args >
      TupleShapeFunctionSet ( Args && ... args )
        : shapeFunctionSetTuple_( std::forward< Args >( args ) ... )
      {
        offset_[ 0 ] = 0;
        Fem::ForLoop< Offsets, 0, sizeof ... ( ShapeFunctionSets ) -1 >::apply( shapeFunctionSetTuple_, offset_ );
      }

      explicit TupleShapeFunctionSet ( const std::tuple< ShapeFunctionSets... > &shapeFunctionSetTuple  )
        : shapeFunctionSetTuple_( shapeFunctionSetTuple )
      {
        offset_[ 0 ] = 0;
        Fem::ForLoop< Offsets, 0, sizeof ... ( ShapeFunctionSets ) -1 >::apply( shapeFunctionSetTuple_, offset_ );
      }

      int order () const { return order( std::index_sequence_for< ShapeFunctionSets ... >() ); }

      std::size_t size () const { return size( std::index_sequence_for< ShapeFunctionSets ... >() ); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        Fem::ForLoop< EvaluateEach, 0, sizeof ... ( ShapeFunctionSets ) -1 >::apply( shapeFunctionSetTuple_, offset_, x, functor );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        Fem::ForLoop< JacobianEach, 0, sizeof ... ( ShapeFunctionSets ) -1 >::apply( shapeFunctionSetTuple_, offset_, x, functor );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        Fem::ForLoop< HessianEach, 0, sizeof ... ( ShapeFunctionSets ) -1 >::apply( shapeFunctionSetTuple_, offset_, x, functor );
      }

      template< std::size_t i >
      const SubShapeFunctionSetType< i > &subShapeFunctionSet ( std::integral_constant< std::size_t, i > = {} ) const
      {
        return std::get< i >( shapeFunctionSetTuple_ );
      }

    protected:
      template< std::size_t ... I >
      int order ( std::index_sequence< I ... > ) const
      {
        return Std::max( std::get< I >( shapeFunctionSetTuple_ ).order() ... );
      }

      template< std::size_t ... I >
      std::size_t size ( std::index_sequence< I ... > ) const
      {
        return Std::sum( std::get< I >( shapeFunctionSetTuple_ ).size() ... );
      }

      template< int >
      static GeometryType makeGeometryType ( GeometryType type )
      {
        return type;
      }

      template< std::size_t ... I >
      static std::tuple< decltype( makeGeometryType< I >( std::declval< GeometryType >() ) ) ... >
      makeGeometryTypeTuple ( GeometryType type, std::index_sequence< I ... > )
      {
        return std::make_tuple( makeGeometryType< I >( type ) ... );
      }

      std::tuple< ShapeFunctionSets... > shapeFunctionSetTuple_;
      Offset offset_;
    };



    // TupleShapeFunctionSet::Offsets
    // ------------------------------

    template< class ... ShapeFunctionSets >
    template< int I >
    struct TupleShapeFunctionSet< ShapeFunctionSets ... >::Offsets
    {
      template< class Tuple >
      static void apply ( const Tuple &tuple, Offset &offset )
      {
        offset[ I + 1 ] = offset[ I ] + std::get< I >( tuple ).size();
      }
    };



    // TupleShapeFunctionSet::FunctorWrapper
    // -------------------------------------

    template< class ... ShapeFunctionSets >
    template< class Functor, class Value, int I >
    struct TupleShapeFunctionSet< ShapeFunctionSets ... >::FunctorWrapper
    {
      static const int rangeOffset = RangeOffsets< ShapeFunctionSets::FunctionSpaceType::dimRange ... >::template offset< I >();

      explicit FunctorWrapper ( const Functor &functor, const Offset &offset )
        : functor_( functor ), offset_( offset ) {}

      template< class Scalar >
      void operator() ( const std::size_t i, const Scalar &subValue )
      {
        Value value( typename FieldTraits< Value >::field_type( 0.0 ) );
        std::copy( subValue.begin(), subValue.end(), value.begin() + rangeOffset );
        functor_( offset_[ I ]  + i, value );
      }

      template< class K >
      void operator () ( const std::size_t i, const Dune::FieldVector< K, 1 > &subValue )
      {
        MakeVectorialExpression< Dune::FieldVector< K, 1 >, Value > value( rangeOffset, subValue[ 0 ] );
        functor_( offset_[ I ]  + i, value );
      }

      template< class K, int n >
      void operator () ( const std::size_t i, const Dune::FieldMatrix< K, 1, n > &subValue )
      {
        MakeVectorialExpression< Dune::FieldMatrix< K, 1, n >, Value > value( rangeOffset, subValue[ 0 ] );
        functor_( offset_[ I ]  + i, value );
      }

      template< class Scalar, class Vectorial >
      void operator() ( const std::size_t i, const MakeVectorialExpression< Scalar, Vectorial > &subValue )
      {
        MakeVectorialExpression< Scalar, Value > value( subValue.component() + rangeOffset, subValue.scalar() );
        functor_( offset_[ I ]  + i, value );
      }

    private:
      Functor functor_;
      const Offset &offset_;
    };



    // TupleShapeFunctionSet::EvaluateEach
    // -----------------------------------

    template< class ... ShapeFunctionSets >
    template< int I >
    struct TupleShapeFunctionSet< ShapeFunctionSets ... >::EvaluateEach
    {
      template< class Tuple, class Point, class Functor >
      static void apply ( const Tuple &tuple, const Offset &offset, const Point &x, Functor functor )
      {
        FunctorWrapper< Functor, RangeType, I > functorWrapper( functor, offset );
        std::get< I >( tuple ).evaluateEach( x, functorWrapper );
      }
    };



    // TupleShapeFunctionSet::JacobianEach
    // -----------------------------------

    template< class ... ShapeFunctionSets >
    template< int I >
    struct TupleShapeFunctionSet< ShapeFunctionSets ... >::JacobianEach
    {
      template< class Tuple, class Point, class Functor >
      static void apply ( const Tuple &tuple, const Offset &offset, const Point &x, Functor functor )
      {
        FunctorWrapper< Functor, JacobianRangeType, I > functorWrapper( functor, offset );
        std::get< I >( tuple ).jacobianEach( x, functorWrapper );
      }
    };



    // TupleShapeFunctionSet::HessianEach
    // ----------------------------------

    template< class ... ShapeFunctionSets >
    template< int I >
    struct TupleShapeFunctionSet< ShapeFunctionSets ... >::HessianEach
    {
      template< class Tuple, class Point, class Functor >
      static void apply ( const Tuple &tuple, const Offset &offset, const Point &x, Functor functor )
      {
        FunctorWrapper< Functor, HessianRangeType, I > functorWrapper( functor, offset );
        std::get< I >( tuple ).hessianEach( x, functorWrapper );
      }
    };


  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_TUPLE_HH
