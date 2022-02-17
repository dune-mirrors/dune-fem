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
      typedef typename Traits :: SolverParameterType          SolverParameterType;

      /** \brief true if a preconditioner type is exported and can be set using bind( op, p ) */
      static const bool preconditioningAvailable = true;

      /** \brief default constructor
       *  \note parameter  SolverParameter object to steer the linear solvers
       */
      InverseOperatorInterface( const SolverParameterType& parameter )
        : parameter_( std::make_shared< SolverParameterType >( parameter ) )
        , verbose_( parameter_->verbose() && Dune::Fem::Parameter::verbose() )

      {
        unbind();
      }

      /** \brief application of operator to compute
       *   @f[
       *     w = op^-1( u )
       *   @f].
       *  \param u parameter right hand side of linear problem
       *  \param w initial guess for linear solver
       */
      virtual void operator() ( const DomainFunctionType& u, RangeFunctionType& w ) const
      {
        opApply( u, w );
      }


      /** \brief application of operator to compute
       *   @f[
       *     w = op^-1( u )
       *   @f].
       *  \param u parameter right hand side of linear problem
       *  \param w initial guess for linear solver
       *
       *  \note Calling the inverse operator for arbitrary discrete functions
       *  a copy to solver compatible discrete function is made.
       */
      template <class DImpl, class RImpl>
      void operator() ( const DiscreteFunctionInterface< DImpl >&u,
                        DiscreteFunctionInterface< RImpl >& w ) const
      {
        opApply( u, w );
      }

      /** \brief store pointer to linear operator
       *  \param op linear operator following the Dune::Fem:Operator interface
       *
       *  \note A dynamic cast to AssembledOperatorType is carried out. For some solvers this is necessary.
       */
      void bind ( const OperatorType &op )
      {
        operator_ = &op;
        assembledOperator_ = dynamic_cast<const AssembledOperatorType*>( &op );
      }

      /** \brief store pointer to linear operator and preconditioner
       *  \param op               linear operator following the Dune::Fem:Operator interface
       *  \param preconditioner   precondition operator
       *
       *  \note A dynamic cast to AssembledOperatorType is carried out. For some solvers this is necessary.
       */
      void bind ( const OperatorType &op, const PreconditionerType &preconditioner )
      {
        bind( op );
        preconditioner_ = &preconditioner;
      }

      /** \brief reset all pointers and internal temporary memory */
      void unbind () { operator_ = nullptr; assembledOperator_ = nullptr; preconditioner_ = nullptr; rhs_.reset(); x_.reset(); }

      /** \brief return number of iterations used in previous call of application operator */
      int iterations () const { return iterations_; }

      /** \brief set number of max linear iterations to be used before an exception is thrown
       *  \param iter  number of max linear iterations
       */
      virtual void setMaxLinearIterations ( const int iter ) {
        parameter_->setMaxIterations( iter );
      }

      /** \brief \copydoc Dune::Fem::InverseOperartorInterface::setMaxLinearIterations */
      virtual void setMaxIterations ( const int iter ) {
        parameter_->setMaxIterations( iter );
      }

      /** \brief set complete set of linear inverse operator parameters
       *  \note newParams  set of new parameters
       */
      void setParameters( const SolverParameterType& newParams)
      {
        std::shared_ptr< SolverParameterType > sharedNewParams = std::make_shared< SolverParameterType > (newParams);
        parameter_.swap( sharedNewParams );
        verbose_ = parameter_->verbose() && Dune::Fem::Parameter::verbose();
      }

      SolverParameterType& parameter () const
      {
        return *parameter_;
      }

      bool verbose() const
      {
        return verbose_;
      }

      //! return accumulated communication time
      double averageCommTime() const
      {
        return -1.;
      }

      //! copy constructor setting defaults
      InverseOperatorInterface(const InverseOperatorInterface &other)
        : parameter_(other.parameter_),
          operator_(nullptr),
          assembledOperator_(nullptr),
          preconditioner_(nullptr),
          rhs_(),
          x_(),
          iterations_(-1),
          rightHandSideCopied_(false)
      {}

    protected:
      // specialization that works with the solvers native storage type
      void opApply( const SolverDiscreteFunctionType& u, SolverDiscreteFunctionType& w ) const
      {
        rightHandSideCopied_ = false;
        iterations_ = asImp().apply( u, w );
      }

      template <class DImpl, class RImpl>
      void opApply( const DiscreteFunctionInterface< DImpl >&u,
                    DiscreteFunctionInterface< RImpl >& w ) const
      {
        if( ! assembledOperator_ )
          DUNE_THROW(Dune::NotImplemented, "InverseOperator::operator() for matrix free operators only makes sense" <<
                                           " for fixed types of domain and range functions to avoid excessive copying!");

        if( ! rhs_ )
        {
          rhs_.reset( new SolverDiscreteFunctionType( "InvOp::rhs", u.space() ) );
        }

        if( ! x_ )
        {
          x_.reset( new SolverDiscreteFunctionType( "InvOp::x", w.space() ) );
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

      std::shared_ptr<SolverParameterType> parameter_;

      const OperatorType*                   operator_ = nullptr;
      const AssembledOperatorType* assembledOperator_ = nullptr;
      const PreconditionerType*       preconditioner_ = nullptr;

      // temporary functions for solver compatibility
      mutable std::unique_ptr< SolverDiscreteFunctionType > rhs_;
      mutable std::unique_ptr< SolverDiscreteFunctionType > x_;

      mutable int  iterations_ = -1 ;
      mutable bool rightHandSideCopied_ = false ;
      mutable bool verbose_;
    };
  } // end namespace Fem
} // end namespace Dune

#endif // DUNE_FEM_SOLVER_INVERSEOPERATORINTERFACE_HH
