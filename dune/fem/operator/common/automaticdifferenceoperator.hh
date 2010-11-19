// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_AUTOMATICDIFFERENCEOPERATOR_HH
#define DUNE_FEM_AUTOMATICDIFFERENCEOPERATOR_HH

#include <limits>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class Operator >
    class AutomaticDifferenceOperator;



    // AutomaticDifferenceLinearOperator
    // ---------------------------------

    template< class Operator >
    class AutomaticDifferenceLinearOperator 
    : public Dune::Fem::Operator< typename Operator::DomainFunctionType, typename Operator::RangeFunctionType >
    {
      typedef Dune::Fem::Operator< typename Operator::DomainFunctionType, typename Operator::RangeFunctionType > BaseType;

      friend class AutomaticDifferenceOperator< Operator >;

    public:
      typedef Operator OperatorType;

      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFieldType RangeFieldType;
      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;

      AutomaticDifferenceLinearOperator ( const std::string &name, const DomainSpaceType &dSpace, const RangeSpaceType &rSpace )
      : name_( name ), 
        op_( 0 ), // initial value for op_ is 'undefined'
        u_( 0 ), // initial value for u_ is 'undefined'
        b_( "AutomaticDifferenceOperator::b_", dSpace ),
        op_u_( "AutomaticDifferenceOperator::op_u_", rSpace ),
        norm_u_( 0 )
      {}

      virtual void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const;

    private:
      void set ( const DomainFunctionType &u, const OperatorType &op, const RangeFieldType &eps );

      const std::string name_;
      const OperatorType *op_;
      const DomainFunctionType* u_;

      mutable DomainFunctionType b_;
      RangeFunctionType op_u_;

      RangeFieldType eps_;
      RangeFieldType norm_u_;
    }; 



    /** \class AutomaticDifferenceOperator
     *  \brief operator wrapper providing a Jacobian through automatic differentiation
     *
     *  \note The Jacobian operator is an on-the-fly operator, i.e., it does
     *        not store a matrix but only implements the application.
     */
    template< class Operator >
    class AutomaticDifferenceOperator 
    : public Dune::Fem::DifferentiableOperator< AutomaticDifferenceLinearOperator< Operator > >
    {
      typedef Dune::Fem::DifferentiableOperator< AutomaticDifferenceLinearOperator< Operator > > BaseType;

    public:
      typedef Operator OperatorType;

      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFieldType RangeFieldType;
      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;

      explicit AutomaticDifferenceOperator ( const Operator &op )
      : operator_(op),
        eps_( Parameter::getValue< RangeFieldType >( "fem.differenceoperator.eps", RangeFieldType( 0 ) ) )
      {}

      AutomaticDifferenceOperator ( const Operator& op, const RangeFieldType &eps )
      : operator_( op ),
        eps_( eps )
      {}

      virtual void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        operator_( arg, dest );
      }

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
      {
        jOp.set( u, operator_, eps_ );
      }
        
    private:
      const Operator &operator_;
      const RangeFieldType eps_;
    };



    // Implementation of AutomaticDifferenceLinearOperator
    // ---------------------------------------------------

    template< class Operator >
    inline void AutomaticDifferenceLinearOperator< Operator >
      ::operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
    {
      assert( op_ && u_ );

      // 'Normal' difference-quotient, i.e.
      // dest = 1/eps (op_(*u +eps*arg) - op_(*u))
      //
      // eps is chosen dynamically, the same way as in 
      // dune-fem/dune/fem/solver/ode/quasi_exact_newton.cpp, see also
      // http://www.freidok.uni-freiburg.de/volltexte/3762/pdf/diss.pdf, page 137
      RangeFieldType eps = eps_;
      if( eps <= RangeFieldType( 0 ) )
      {
        const RangeFieldType machine_eps = std::numeric_limits< RangeFieldType >::epsilon();
        const RangeFieldType norm_p_sq = arg.scalarProductDofs( arg );
        if( norm_p_sq > machine_eps )
          eps = std::sqrt( (RangeFieldType( 1 ) + norm_u_) * machine_eps / norm_p_sq );
        else
          eps = std::sqrt( machine_eps );
      }

      b_.assign( *u_ );
      b_.addScaled( arg, eps );
      (*op_)( b_, dest );
      dest -= op_u_;
      dest *= RangeFieldType( 1 ) / eps;
    }


    template< class Operator >
    inline void AutomaticDifferenceLinearOperator< Operator >
      ::set ( const DomainFunctionType &u, const OperatorType &op, const RangeFieldType &eps )
    { 
      u_ = &u; 
      op_ = &op;
      (*op_)( *u_, op_u_ );

      eps_ = eps;
      if( eps_ <= RangeFieldType( 0 ) )
        norm_u_ = std::sqrt( u_->scalarProductDofs( *u_ ) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_AUTOMATICDIFFERENCEOPERATOR_HH
