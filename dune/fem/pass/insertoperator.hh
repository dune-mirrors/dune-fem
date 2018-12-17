#ifndef DUNE_FEM_PASS_INSERTOPERATOR_HH
#define DUNE_FEM_PASS_INSERTOPERATOR_HH

#if HAVE_DUNE_FEM_DG
#error "Outdated header, #include <dune/fem-dg/pass/insertoperator.hh> instead!"
#endif

#include <memory>
#include <string>
#include <tuple>
#include <type_traits>

#include <dune/common/tupleutility.hh>

#include <dune/fem/common/memory.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/pass/common/pass.hh>

namespace Dune
{

  namespace Fem
  {

#ifndef DOXYGEN

    namespace __InsertOperatorPass
    {

      // isOperator
      // ----------

      template< class DomainFunction, class RangeFunction >
      std::true_type __isOperator ( const Dune::Fem::Operator< DomainFunction, RangeFunction > & );

      std::false_type __isOperator ( ... );

      template< class Operator >
      using isOperator = decltype( __isOperator( std::declval< Operator >() ) );



      template< class Operator >
      struct DiscreteModel
      {
        struct Traits
        {
          typedef typename Operator::RangeFunctionType DestinationType;
        };
      };

    } // namespace __InsertOperatorPass

#endif // #ifndef DOXYGEN



    // InsertOperatorPass
    // ------------------

    /** \brief include a Dune::Fem::Operator into a pass
     *
     *  \note This pass stores a reference to the operator passed to the
     *  constructor. The operator will always be called on the result of
     *  PreviousPass, so the position of this pass in your pass tree is of
     *  importance!
     *
     *  \tparam  Operator  a Dune::Fem::Operator
     *  \tparam  PreviousPass  type of previous pass
     *  \tparam  id  unique id for this pass
     */
    template< class Operator, class PreviousPass, int id >
    class InsertOperatorPass final
      : public Dune::Fem::Pass< __InsertOperatorPass::DiscreteModel< Operator >, PreviousPass, id >
    {
      typedef Dune::Fem::Pass< __InsertOperatorPass::DiscreteModel< Operator >, PreviousPass, id > BaseType;

      static_assert( __InsertOperatorPass::isOperator< Operator >::value, "Template parameter Operator is not derived from Dune::Fem::Operator" );

    public:
      /** \brief export operator type */
      typedef Operator OperatorType;

      /** \copydoc Dune::Fem::Pass::PreviousPassType */
      typedef typename BaseType::PreviousPassType PreviousPassType;

      /** \copydoc Dune::Fem::Pass::TotalArgumentType */
      typedef typename BaseType::TotalArgumentType TotalArgumentType;
      /** \copydoc Dune::Fem::Pass::DestinationType */
      typedef typename BaseType::DestinationType DestinationType;

      /** \copydoc Dune::Fem::Pass::passNumber */
      using BaseType::passNumber;

    protected:
      using BaseType::deleteHandler_;
      using BaseType::destination_;

    public:
      /** \name Construction
       *  \{
       */

      /** \brief constructor
       *
       *  \param[in]  space  range discrete function space
       *  \param[in]  op  an operator instance
       *  \param[in]  pass  previous pass
       */
      InsertOperatorPass ( const typename OperatorType::RangeFunctionType::DiscreteFunctionSpaceType &space,
                           const OperatorType &op, PreviousPassType &pass )
        : BaseType( pass ),
          space_( referenceToSharedPtr( space ) ),
          operator_( op )
      {}

      /** \} */

#ifndef DOXYGEN

      ~InsertOperatorPass () {}

#endif // #ifndef DOXYGEN

      /** \name Public member methods
       *  \{
       */

      /** \copydoc Dune::Fem::Pass::allocateLocalMemory */
      void allocateLocalMemory ()
      {
        if( destination_ )
          return;
        destination_ = new DestinationType( "stabilization_pass_" + std::to_string( passNumber() ), space() );
        deleteHandler_ = &(BaseType::DeleteHandlerType::instance());
      }

      /** \copydoc Dune::Fem::Pass::printTexInfo */
      void printTexInfo ( std::ostream & ) const
      {
        DUNE_THROW( NotImplemented, "Method printTexInfo() not implemented" );
      }

      /** \copydoc Dune::Fem::Pass::compute */
      void compute ( const TotalArgumentType &arg, DestinationType &dest ) const
      {
        apply( arg, dest );
      }

      /** \} */

      /** \name Non-interface methods
       *  \{
       */

      /** \brief return discrete function space */
      const typename OperatorType::RangeFunctionType::DiscreteFunctionSpaceType &space () const
      {
        return *space_;
      }

      /** \} */

    protected:
      void apply ( const TotalArgumentType &u, DestinationType &w ) const
      {
        // get argument
        typedef std::integral_constant< int, PreviousPassType::passId > Key;
        const std::size_t position = Dune::FirstTypeIndex< typename PreviousPassType::PassIds, Key >::type::value;
        const typename OperatorType::DomainFunctionType &v = *std::get< position >( u );

        // apply operator
        operator_( v, w );
      }

    private:
      std::shared_ptr< const typename OperatorType::RangeFunctionType::DiscreteFunctionSpaceType > space_;
      const OperatorType &operator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PASS_INSERTOPERATOR_HH
