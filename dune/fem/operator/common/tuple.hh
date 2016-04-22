#ifndef DUNE_FEM_OPERATOR_COMMON_TUPLE_HH
#define DUNE_FEM_OPERATOR_COMMON_TUPLE_HH

#include <tuple>
#include <utility>

#include <dune/fem/common/tupleforeach.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/operator/common/operator.hh>


namespace Dune
{

  namespace Fem
  {

    // forward declaration
    template< class ... Operators >
    class TupleOperator;


    namespace __TupleOperatorImp
    {

      // RowTraits
      // ---------

      template< class ... Operators >
      struct RowTraits
      {
        typedef TupleDiscreteFunction< typename Operators::DomainFunctionType ... > DomainFunctionType;
        static_assert( Std::are_all_same< typename Operators::RangeFunctionType ... >::value,
                       "TupleOperator< RowTraits > needs a common RangeFunction Type." );
        typedef typename std::tuple_element< 0, std::tuple< Operators ... > >::type::RangeFunctionType RangeFunctionType;
      };

      // ColTraits
      // ---------

      template< class ... Operators >
      struct ColTraits
      {
        static_assert( Std::are_all_same< typename Operators::DomainFunctionType ... >::value,
                       "TupleOperator< ColTraits > needs a common DomainFunction Type." );
        typedef typename std::tuple_element< 0, std::tuple< Operators ... > >::type::DomainFunctionType DomainFunctionType;
        typedef TupleDiscreteFunction< typename Operators::RangeFunctionType ... > RangeFunctionType;
      };

      template< class ... Operators >
      struct Traits : public RowTraits< Operators ... > {};

      template< class ... Operators >
      struct Traits< std::tuple< Operators > ... > : public ColTraits< Operators ... > {};

    }


    // TupleOperator
    // -------------

    template< class ... Operators >
    class TupleOperator
      : public Operator< typename __TupleOperatorImp::Traits< Operators ... >::DomainFunctionType,
                         typename __TupleOperatorImp::Traits< Operators ... >::RangeFunctionType >,
        public std::tuple< Operators ... >
    {
      typedef TupleOperator< Operators ... > ThisType;
      typedef typename __TupleOperatorImp::Traits< Operators ... > Traits;
      typedef Operator< typename Traits::DomainFunctionType, typename Traits::RangeFunctionType > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      template< class ... Args >
      TupleOperator ( std::tuple< Args ... > && args )
        : BaseType( std::forward( args ) )
      {}

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        dest.clear();
        RangeFunctionType tmp( dest );
        auto functor = [ &arg, &dest, &tmp ] ( auto &op, auto I )
        {
          op( std::get< I >( arg ), tmp );
          dest += tmp;
        };

        for_each( *this, functor );
      }
    };


    template< class ... Operators >
    class TupleOperator< std::tuple< Operators > ... >
      : public Operator< typename __TupleOperatorImp::Traits< std::tuple< Operators > ... >::DomainFunctionType,
                         typename __TupleOperatorImp::Traits< std::tuple< Operators > ... >::RangeFunctionType >,
        public std::tuple< Operators ... >
    {
      typedef TupleOperator< Operators ... > ThisType;
      typedef typename __TupleOperatorImp::Traits< std::tuple< Operators > ... > Traits;
      typedef Operator< typename Traits::DomainFunctionType, typename Traits::RangeFunctionType > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      template< class ... Args >
      TupleOperator ( std::tuple< Args ... > && args )
        : BaseType( std::forward( args ) )
      {}

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        dest.clear();
        auto functor = [ &arg, &dest ] ( auto &op, auto I ) { op( arg, std::get< I >( dest ) ); };

        for_each( *this, functor );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_TUPLE_HH
