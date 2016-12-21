#ifndef DUNE_FEMPY_PY_OPERATOR_HH
#define DUNE_FEMPY_PY_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/corepy/pybind11/pybind11.h>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // GeneralGridFunction
      // -------------------

      template< class DiscreteFunction >
      using GeneralGridFunction = VirtualizedGridFunction< typename DiscreteFunction::GridPartType, typename DiscreteFunction::RangeType >;



      // registerGeneralOperatorCall
      // ---------------------------

      template< class Operator, class Holder, class Alias, decltype( std::declval< const Operator & >()( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::RangeFunctionType & >() ), 0 ) = 0 >
      void registerGeneralOperatorCall ( pybind11::class_< Operator, Holder, Alias > &cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "__call__", [] ( Operator &op, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::RangeFunctionType &w ) { op( u, w ); }, "u"_a, "w"_a );
      }

      template< class Operator, class Holder, class Alias >
      void registerOperatorGeneralCall ( pybind11::class_< Operator, Holder, Alias > &cls, PriorityTag< 0 > )
      {}



      // registerGeneralOperatorJacobian
      // -------------------------------

      template< class Operator, class Holder, class Alias, decltype( std::declval< const Operator & >().jacobian( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::JacobianOperatorType >() ), 0 ) = 0 >
      void registerGeneralOperatorJacobian ( pybind11::class_< Operator, Holder, Alias > &cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "jacobian", [] ( Operator &op, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::JacobianRangeType &jOp ) { op.jacobian( u, jOp ); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class Holder, class Alias >
      void registerGeneralOperatorJacobian ( pybind11::class_< Operator, Holder, Alias > &cls, PriorityTag< 0 > )
      {}



      // registerOperatorJacobian
      // ------------------------

      template< class Operator, class Holder, class Alias, decltype( std::declval< const Operator & >().jacobian( std::declval< const typename Operator::DomainFunctionType & >(), std::declval< typename Operator::JacobianOperatorType >() ), 0 ) = 0 >
      void registerOperatorJacobian ( pybind11::class_< Operator, Holder, Alias > &cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "jacobian", [] ( Operator &op, const typename Operator::DomainFunctionType &u, typename Operator::JacobianRangeType &jOp ) { op.jacobian( u, jOp ); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class Holder, class Alias >
      void registerOperatorJacobian ( pybind11::class_< Operator, Holder, Alias > &cls, PriorityTag< 0 > )
      {}



      // registerOperator
      // ----------------

      template< class Operator, class Cls >
      void registerOperator ( pybind11::module module, Cls &cls )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename Operator::RangeFunctionType RangeFunction;

        using pybind11::operator""_a;

        cls.def( "__call__", [] ( Operator &op, const DomainFunction &u, RangeFunction &w ) { op( u, w ); }, "u"_a, "w"_a );
        registerGeneralOperatorCall( cls, PriorityTag< 42 >() );

        registerOperatorJacobian( cls, PriorityTag< 42 >() );
        registerGeneralOperatorJacobian( cls, PriorityTag< 42 >() );
      }

    } // namespace detail



    // registerOperator
    // ----------------

    template< class Operator, class Holder, class Alias >
    void registerOperator ( pybind11::module module, pybind11::class_< Operator, Holder, Alias > &cls )
    {
      detail::registerOperator< Operator >( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_OPERATOR_HH
