#ifndef DUNE_FEM_OPERATOR_COMMON_TUPLE_HH
#define DUNE_FEM_OPERATOR_COMMON_TUPLE_HH

#include <tuple>
#include <utility>

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

    }


    // TupleOperator
    // -------------

    template< class ... Operators >
    class TupleOperator
      : public Operator< typename __TupleOperatorImp::RowTraits< Operators ... >::DomainFunctionType,
                         typename __TupleOperatorImp::RowTraits< Operators ... >::RangeFunctionType >,
        public std::tuple< Operators ... >
    {
      typedef TupleOperator< Operators ... > ThisType;
      // Note we have to check the rows of each column for consistency
      typedef typename __TupleOperatorImp::RowTraits< Operators ... > Traits;
      typedef Operator< typename Traits::DomainFunctionType, typename Traits::RangeFunctionType > BaseType;
      typedef std::tuple< Operators ... > TupleType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      template< class ... Args >
      TupleOperator ( Args&& ... args )
        : TupleType( std::forward< Args >( args ) ... )
      {}

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        dest.clear();
        RangeFunctionType tmp( dest );
        apply( arg, dest, tmp, std::integral_constant< std::size_t, 0 >() );
      }

    protected:
      template< std::size_t I >
      void apply ( const DomainFunctionType &arg, RangeFunctionType &dest,
                   RangeFunctionType &tmp, std::integral_constant< std::size_t, I > ) const
      {
        std::get< I >( *this )( std::get< I >( arg ), tmp );
        dest += tmp;
        apply( arg, dest, tmp, std::integral_constant< std::size_t, I + 1 >() );
      }

      void apply ( const DomainFunctionType &arg, RangeFunctionType &dest,
                   RangeFunctionType &tmp, std::integral_constant< std::size_t, sizeof ... (Operators ) > ) const
      {}
    };



    // RowTupleOperator
    // ----------------

    template< class ... Operators >
    class RowTupleOperator
      : public Operator< typename __TupleOperatorImp::ColTraits< Operators ... >::DomainFunctionType,
                         typename __TupleOperatorImp::ColTraits< Operators ... >::RangeFunctionType >,
        public std::tuple< Operators ... >
    {
      typedef RowTupleOperator< Operators ... > ThisType;
      // Note we have to check the columns of each row for consistency
      typedef typename __TupleOperatorImp::ColTraits< Operators ... > Traits;
      typedef Operator< typename Traits::DomainFunctionType, typename Traits::RangeFunctionType > BaseType;
      typedef std::tuple< Operators ... > TupleType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      template< class ... Args >
      RowTupleOperator ( Args&& ... args )
        : TupleType( std::forward< Args >( args ) ... )
      {}

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        dest.clear();
        apply( arg, dest, std::integral_constant< std::size_t, 0 >() );
      }

    protected:
      template< std::size_t I >
      void apply ( const DomainFunctionType &arg, RangeFunctionType &dest, std::integral_constant< std::size_t, I > ) const
      {
        std::get< I >( *this )( arg,  std::get< I >( dest ) );
        apply( arg, dest, std::integral_constant< std::size_t, I + 1 >() );
      }

      void apply ( const DomainFunctionType &arg, RangeFunctionType &dest, std::integral_constant< std::size_t, sizeof ... (Operators ) > ) const
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_TUPLE_HH
