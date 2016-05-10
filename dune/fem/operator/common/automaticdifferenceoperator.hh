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

    template< class DomainFunction, class RangeFunction, class LinearOperator >
    class AutomaticDifferenceOperator;



    // AutomaticDifferenceLinearOperator
    // ---------------------------------

    template< class DomainFunction, class RangeFunction = DomainFunction >
    class AutomaticDifferenceLinearOperator
    : public Dune::Fem::Operator< DomainFunction, RangeFunction >
    {
      typedef AutomaticDifferenceLinearOperator< DomainFunction, RangeFunction > ThisType;
      typedef Dune::Fem::Operator< DomainFunction, RangeFunction > BaseType;

      friend class AutomaticDifferenceOperator< DomainFunction, RangeFunction, ThisType >;

    public:
      typedef Dune::Fem::Operator< DomainFunction, RangeFunction > OperatorType;

      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFieldType RangeFieldType;
      typedef typename BaseType::DomainFieldType DomainFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

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

    protected:
      void set ( const DomainFunctionType &u, const OperatorType &op, const RealType &eps );
      const std::string name_;
      const OperatorType *op_;
      const DomainFunctionType *u_;

      mutable DomainFunctionType b_;
      RangeFunctionType op_u_;

      RangeFieldType eps_;
      RangeFieldType norm_u_;
    };



    // AutomaticDifferenceOperator
    // ---------------------------

    /** \class AutomaticDifferenceOperator
     *  \brief operator providing a Jacobian through automatic differentiation
     *
     *  \note The Jacobian operator is an on-the-fly operator, i.e., it does
     *        not store a matrix but only implements the application.
     */
    template< class DomainFunction, class RangeFunction = DomainFunction,
              class LinearOperator = AutomaticDifferenceLinearOperator< DomainFunction, RangeFunction > >
    class AutomaticDifferenceOperator
    : public Dune::Fem::DifferentiableOperator< LinearOperator >
    {
      typedef Dune::Fem::DifferentiableOperator< LinearOperator > BaseType;

    public:
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFieldType RangeFieldType;
      typedef typename BaseType::DomainFieldType DomainFieldType;
      typedef typename Dune::FieldTraits< RangeFieldType >::real_type RealType;

      typedef typename BaseType::JacobianOperatorType JacobianOperatorType;

      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;

      AutomaticDifferenceOperator ( const ParameterReader &parameter = Parameter::container() )
      : eps_( parameter.getValue< RangeFieldType >( "fem.differenceoperator.eps", RangeFieldType( 0 ) ) )
      {}

      explicit AutomaticDifferenceOperator ( const RangeFieldType &eps )
      : eps_( eps )
      {}

      virtual void jacobian ( const DomainFunctionType &u, JacobianOperatorType &jOp ) const
      {
        jOp.set( u, *this, eps_ );
      }

    private:
      const RealType eps_;
    };



    // Implementation of AutomaticDifferenceLinearOperator
    // ---------------------------------------------------

    template< class DomainFunction, class RangeFunction >
    inline void AutomaticDifferenceLinearOperator< DomainFunction, RangeFunction >
      ::operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
    {
      assert( op_ && u_ );

      // 'Normal' difference-quotient, i.e.
      // dest = 1/eps (op_(*u +eps*arg) - op_(*u))
      //
      // eps is chosen dynamically, the same way as in
      // dune-fem/dune/fem/solver/ode/quasi_exact_newton.cpp, see also
      // http://www.freidok.uni-freiburg.de/volltexte/3762/pdf/diss.pdf, page 137
      RealType eps = eps_;
      if( eps <= RealType( 0. ) )
      {
        const RealType machine_eps = std::numeric_limits< RealType >::epsilon();
        const RealType norm_p_sq = arg.normSquaredDofs( );
        if( norm_p_sq > machine_eps )
          eps = std::sqrt( (RealType( 1 ) + norm_u_) * machine_eps / norm_p_sq );
        else
          eps = std::sqrt( machine_eps );
      }

      b_.assign( *u_ );
      b_.axpy( eps, arg );
      (*op_)( b_, dest );
      dest -= op_u_;
      dest *= RealType( 1 ) / eps;
    }


    template< class DomainFunction, class RangeFunction >
    inline void AutomaticDifferenceLinearOperator< DomainFunction, RangeFunction >
      ::set ( const DomainFunctionType &u, const OperatorType &op, const RealType &eps )
    {
      u_ = &u;
      op_ = &op;
      (*op_)( *u_, op_u_ );

      eps_ = eps;
      if( eps_ <= RealType( 0 ) )
        norm_u_ = std::sqrt( u_->scalarProductDofs( *u_ ) );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_AUTOMATICDIFFERENCEOPERATOR_HH
