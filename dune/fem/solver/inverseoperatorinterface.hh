#ifndef DUNE_FEM_SOLVER_INVERSEOPERATORINTERFACE_HH
#define DUNE_FEM_SOLVER_INVERSEOPERATORINTERFACE_HH

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/parameter.hh>
#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune {
  namespace Fem {

    template <class Traits>
    class InverseOperatorInterface :
      public Dune::Fem::Operator< typename Traits::DiscreteFunctionType, typename Traits::DiscreteFunctionType >,
      public BartonNackmanInterface< InverseOperatorInterface< Traits >, typename Traits::InverseOperatorType >
    {
    protected:
      typedef typename Traits::OperatorType         BaseType;
      typedef BartonNackmanInterface< InverseOperatorInterface< Traits >, typename Traits::InverseOperatorType > Base2Type;

      using Base2Type :: asImp;

      typedef typename Traits::InverseOperatorType  InverseOperatorType;
    public:
      typedef typename BaseType :: DomainFunctionType DomainFunctionType;
      typedef typename BaseType :: RangeFunctionType  RangeFunctionType;

      typedef typename Traits :: SolverDiscreteFunctionType   SolverDiscreteFunctionType;
      typedef typename Traits :: OperatorType                 OperatorType;
      typedef typename Traits :: AssembledOperatorType        AssembledOperatorType;
      typedef typename Traits :: PreconditionerType           PreconditionerType;

      InverseOperatorInterface( const SolverParameter& parameter = SolverParameter(Parameter::container() ) )
        : parameter_( std::make_shared< SolverParameter >( parameter ) )
      {}

      //! application of operator, i.e. solution of inverse operator with given right hand side and initial guess
      // TODO: improve docu
      virtual void operator() ( const SolverDiscreteFunctionType& u, SolverDiscreteFunctionType& w ) const
      {
        rightHandSideCopied_ = false;
        iterations_ = asImp().apply( u, w );
      }

      //! application of operator, here solution and right hand side can be of any discrete function type
      // TODO: improve docu
      template <class DImpl, class RImpl>
      void operator() ( const DiscreteFunctionInterface< DImpl >&u,
                        DiscreteFunctionInterface< RImpl >& w ) const
      {
        if( ! assembledOperator_ )
          DUNE_THROW(Dune::NotImplemented, "InverseOperator::operator() for matrix free operators only makes sense" <<
                                           " for fixed types of domain and range functions to avoid excessive copying!");

        if( ! rhs_ )
        {
          rhs_.reset( new DomainFunctionType( "InvOp::rhs", u.space() ) );
        }

        if( ! x_ )
        {
          x_.reset( new RangeFunctionType( "InvOp::x", w.space() ) );
        }

        // copy right hand side
        rhs_->assign( u );
        rightHandSideCopied_ = true;

        // copy initial guess
        x_->assign( w );

        iterations_ = asImp().apply( *rhs_, *x_ );

        // store result in destination
        w.assign( *x_ );
        rightHandSideCopied_ = false;
      }

      void bind ( const OperatorType &op )
      {
        operator_ = &op;
        assembledOperator_ = dynamic_cast<const AssembledOperatorType*>( &op );
      }

      void bind ( const OperatorType &op, const PreconditionerType &preconditioner )
      {
        bind( op );
        preconditioner_ = &preconditioner;
      }

      void unbind () { operator_ = nullptr; assembledOperator_ = nullptr; preconditioner_ = nullptr; rhs_.reset(); x_.reset(); }

      int iterations () const { return iterations_; }

      virtual void setMaxLinearIterations ( const int iter ) {
        parameter_->setMaxIterations( iter );
      }

      virtual void setMaxIterations ( const int iter ) {
        parameter_->setMaxIterations( iter );
      }

      void setParameters( const SolverParameter& newParams)
      {
        std::shared_ptr< SolverParameter > sharedNewParams = std::make_shared< SolverParameter > (newParams);
        parameter_.swap( sharedNewParams );
      }

      SolverParameter& parameter () const
      {
        return *parameter_;
      }

      //! return accumulated communication time
      double averageCommTime() const
      {
        return -1;
      }

    protected:
      std::shared_ptr<SolverParameter> parameter_;

      const OperatorType*                   operator_ = nullptr;
      const AssembledOperatorType* assembledOperator_ = nullptr;
      const PreconditionerType*       preconditioner_ = nullptr;

      // temporary functions for solver compatibility
      mutable std::unique_ptr< DomainFunctionType > rhs_;
      mutable std::unique_ptr< RangeFunctionType  > x_;

      mutable int iterations_ = -1 ;

      mutable bool rightHandSideCopied_ = false ;
    };
  } // end namespace Fem
} // end namespace Dune

#endif // DUNE_FEM_SOLVER_INVERSEOPERATORINTERFACE_HH
