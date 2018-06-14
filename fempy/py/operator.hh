#ifndef DUNE_FEMPY_PY_OPERATOR_HH
#define DUNE_FEMPY_PY_OPERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/common/typeutilities.hh>

#include <dune/fem/operator/common/differentiableoperator.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/pybind11/pybind11.hh>

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

      template< class Operator, class... options, decltype( std::declval< const Operator & >()( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::RangeFunctionType & >() ), 0 ) = 0 >
      inline static void registerGeneralOperatorCall ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "__call__", [] ( Operator &self, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::RangeFunctionType &w ) { self( u, w ); }, "u"_a, "w"_a );
      }

      template< class Operator, class... options >
      inline static void registerGeneralOperatorCall ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}



      // registerGeneralOperatorJacobian
      // -------------------------------

      template< class Operator, class... options, decltype( std::declval< const Operator & >().jacobian( std::declval< const GeneralGridFunction< typename Operator::DomainFunctionType > & >(), std::declval< typename Operator::JacobianOperatorType >() ), 0 ) = 0 >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;
        cls.def( "jacobian", [] ( Operator &self, const GeneralGridFunction< typename Operator::DomainFunctionType > &u, typename Operator::JacobianRangeType &jOp ) { self.jacobian( u, jOp ); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static void registerGeneralOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}



      // registerOperatorJacobian
      // ------------------------

      template< class Operator, class... options, decltype( std::declval< const Operator & >().jacobian( std::declval< const typename Operator::DomainFunctionType & >(), std::declval< typename Operator::JacobianOperatorType >() ), 0 ) = 0 >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 1 > )
      {
        using pybind11::operator""_a;

        cls.def( "jacobian", [] ( Operator &self, const typename Operator::DomainFunctionType &u, typename Operator::JacobianRangeType &jOp ) { self.jacobian( u, jOp ); }, "u"_a, "jOp"_a );
      }

      template< class Operator, class... options >
      inline static void registerOperatorJacobian ( pybind11::class_< Operator, options... > cls, PriorityTag< 0 > )
      {}



      // registerOperator
      // ----------------

      template< class Operator, class... options >
      inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
      {
        typedef typename Operator::DomainFunctionType DomainFunction;
        typedef typename Operator::RangeFunctionType RangeFunction;

        using pybind11::operator""_a;

        cls.def( "__call__", [] ( Operator &self, const DomainFunction &u, RangeFunction &w ) { self( u, w ); }, "u"_a, "w"_a );
        registerGeneralOperatorCall( cls, PriorityTag< 42 >() );

        registerOperatorJacobian( cls, PriorityTag< 42 >() );
        registerGeneralOperatorJacobian( cls, PriorityTag< 42 >() );
      }

    } // namespace detail



    // registerOperator
    // ----------------

    template< class Operator, class... options >
    inline static void registerOperator ( pybind11::module module, pybind11::class_< Operator, options... > cls )
    {
      detail::registerOperator< Operator >( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_OPERATOR_HH
