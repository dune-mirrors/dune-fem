#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH

#include <cassert>
#include <cmath>

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <utility>

#include <dune/fem/solver/parameter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class EisenstatWalkerStrategy
     *
     * \brief Adaptive tolerance selection for linear solver.
     *
     * \note Prevents oversolving linear systems far away from solution.
     *       Source: "Globally Convergent Inexact Newton Methods", Stanley C. Eisenstat and Homer F. Walker, https://doi.org/10.1137/0804022.
    */
    class EisenstatWalkerStrategy
    {
    protected:
      const double etaMax_ = 0.99;
      const double gamma_ = 0.1;
      mutable double previousEta_ = -1.0;
      mutable double previousResidual_ = -1.0;
      mutable double newtonTolerance_;

    public:
      /** constructor
       *  \param[in]  newtonTolerance      the absolute tolerance of the Newton method
      */
      EisenstatWalkerStrategy(const double newtonTolerance) : newtonTolerance_(newtonTolerance) {}
      double nextLinearTolerance(const double currentResidual) const
      {
        double eta = etaMax_;
        // First call previousEta_ is negative
        if (previousEta_ >= 0.0)
        {
          const double etaA = gamma_ * currentResidual * currentResidual / (previousResidual_ * previousResidual_);
          const double indicator = gamma_ * previousEta_ * previousEta_;
          const double etaC = indicator < 0.1 ? std::min(etaA, etaMax_) : std::min(etaMax_, std::max(etaA, indicator));
          eta = std::min(etaMax_, std::max(etaC, 0.5 * newtonTolerance_ / currentResidual));
        }
        previousResidual_ = currentResidual;
        previousEta_ = eta;
        return eta;
      }
      void setTolerance(const double newtonTolerance) { newtonTolerance_ = newtonTolerance; }
    };

    // NewtonParameter
    // ---------------

    template <class SolverParam = SolverParameter>
    struct NewtonParameter
      : public Dune::Fem::LocalParameter< NewtonParameter<SolverParam>, NewtonParameter<SolverParam> >
    {
    protected:

      std::shared_ptr<SolverParam> baseParam_;
      // key prefix, default is fem.solver.newton. (can be overloaded by user)
      const std::string keyPrefix_;

      ParameterReader parameter_;

    public:
      NewtonParameter( const SolverParam& baseParameter, const std::string keyPrefix = "fem.solver.newton." )
        : baseParam_( static_cast< SolverParam* > (baseParameter.clone()) ),
          keyPrefix_( keyPrefix ),
          parameter_( baseParameter.parameter() )
      {}

      template <class Parameter, std::enable_if_t<!std::is_base_of<SolverParam,Parameter>::value && !std::is_same<Parameter,ParameterReader>::value,int> i=0>
      NewtonParameter( const Parameter& solverParameter, const std::string keyPrefix = "fem.solver.newton." )
        : baseParam_( new SolverParam(solverParameter) ),
          keyPrefix_( keyPrefix ),
          parameter_( solverParameter.parameter() )
      {}

      template <class ParamReader, std::enable_if_t<!std::is_same<ParamReader,SolverParam>::value && std::is_same<ParamReader,ParameterReader>::value,int> i=0>
      NewtonParameter( const ParamReader &parameter, const std::string keyPrefix = "fem.solver.newton." )
        : baseParam_( std::make_shared<SolverParam>( keyPrefix + "linear.", parameter) ),
          keyPrefix_( keyPrefix),
          parameter_( parameter )
      {}

      const ParameterReader &parameter () const { return parameter_; }
      const SolverParam& solverParameter () const { return *baseParam_; }
      const SolverParam& linear () const { return *baseParam_; }

      virtual void reset ()
      {
        baseParam_->reset();
        tolerance_ = -1;
        verbose_ = -1;
        maxIterations_ = -1;
        maxLinearIterations_ = -1;
        maxLineSearchIterations_ = -1;
      }

      //These methods affect the nonlinear solver
      virtual double tolerance () const
      {
        if(tolerance_ < 0)
          tolerance_ =  parameter_.getValue< double >( keyPrefix_ + "tolerance", 1e-6 );
        return tolerance_;
      }

      virtual void setTolerance ( const double tol )
      {
        assert( tol > 0 );
        tolerance_ = tol;
      }

      virtual bool verbose () const
      {
        if(verbose_ < 0)
        {
          // the following causes problems with different default values
          // used if baseParam_->keyPrefix is not default but the default is
          // also used in the program
          // const bool v = baseParam_? baseParam_->verbose() : false;
          const bool v = false;
          verbose_ = parameter_.getValue< bool >(keyPrefix_ +  "verbose", v ) ? 1 : 0 ;
        }
        return verbose_;
      }

      virtual void setVerbose( bool verb)
      {
        verbose_ = verb ? 1 : 0;
      }

      virtual int maxIterations () const
      {
        if(maxIterations_ < 0)
          maxIterations_ =  parameter_.getValue< int >( keyPrefix_ + "maxiterations", std::numeric_limits< int >::max() );
        return maxIterations_;
      }

      virtual void setMaxIterations ( const int maxIter )
      {
        assert(maxIter >= 0);
        maxIterations_ = maxIter;
      }

      //Maximum Linear Iterations in total
      //!= max iterations of each linear solve
      virtual int maxLinearIterations () const
      {
        if(maxLinearIterations_ < 0)
          maxLinearIterations_ = linear().maxIterations();
        return maxLinearIterations_;
      }

      virtual void setMaxLinearIterations ( const int maxLinearIter )
      {
        assert(maxLinearIter >=0);
        maxLinearIterations_ = maxLinearIter;
      }

      virtual int maxLineSearchIterations () const
      {
        if(maxLineSearchIterations_ < 0)
          maxLineSearchIterations_ = parameter_.getValue< int >( keyPrefix_ + "maxlinesearchiterations", std::numeric_limits< int >::max() );
        return maxLineSearchIterations_;
      }

      virtual void setMaxLineSearchIterations ( const int maxLineSearchIter )
      {
        assert( maxLineSearchIter >= 0);
        maxLineSearchIterations_ = maxLineSearchIter;
      }

      enum class LineSearchMethod {
          none   = 0,
          simple = 1
        };

      virtual LineSearchMethod lineSearch () const
      {
        const std::string lineSearchMethods[] = { "none", "simple" };
        return static_cast< LineSearchMethod>( parameter_.getEnum( keyPrefix_ + "lineSearch", lineSearchMethods, 0 ) );
      }

      virtual void setLineSearch ( const LineSearchMethod method )
      {
        const std::string lineSearchMethods[] = { "none", "simple" };
        Parameter::append( keyPrefix_ + "lineSearch", lineSearchMethods[int(method)], true );
      }

      enum class LinearToleranceStrategy {
        none = 0,
        eisenstatwalker = 1
      };

      virtual LinearToleranceStrategy linearToleranceStrategy () const
      {
        const std::string linearToleranceStrategy[] = { "none", "eisenstatwalker" };
        return static_cast< LinearToleranceStrategy>( parameter_.getEnum( keyPrefix_ + "linear.tolerance.strategy", linearToleranceStrategy, 0 ) );
      }

      virtual void setLinearToleranceStrategy ( const LinearToleranceStrategy strategy )
      {
        const std::string linearToleranceStrategy[] = { "none", "eisenstatwalker" };
        Parameter::append( keyPrefix_ + "linear.tolerance.strategy", linearToleranceStrategy[int(strategy)], true );
      }

    private:
      mutable double tolerance_ = -1;
      mutable int verbose_ = -1;
      mutable int maxIterations_ = -1;
      mutable int maxLinearIterations_ = -1;
      mutable int maxLineSearchIterations_ = -1;
    };



    // NewtonFailure
    // -------------

    enum class NewtonFailure
    // : public int
    {
      Success = 0,
      InvalidResidual = 1,
      IterationsExceeded = 2,
      LinearIterationsExceeded = 3,
      LineSearchFailed = 4,
      TooManyIterations = 5,
      TooManyLinearIterations = 6,
      LinearSolverFailed = 7
    };



    // NewtonInverseOperator
    // ---------------------

    /** \class NewtonInverseOperator
     *  \brief inverse operator based on a newton scheme
     *
     *  \tparam  Op      operator to invert (must be a DifferentiableOperator)
     *  \tparam  LInvOp  linear inverse operator
     *
     *  \note Verbosity of the NewtonInverseOperator is controlled via the
     *        paramter <b>fem.solver.newton.verbose</b>; it defaults to
     *        <b>fem.solver.verbose</b>.
     */
    template< class JacobianOperator, class LInvOp >
    class NewtonInverseOperator
    : public Operator< typename JacobianOperator::RangeFunctionType, typename JacobianOperator::DomainFunctionType >
    {
      typedef NewtonInverseOperator< JacobianOperator, LInvOp > ThisType;
      typedef Operator< typename JacobianOperator::RangeFunctionType, typename JacobianOperator::DomainFunctionType > BaseType;

      //! selects LInvOp::preconditioningAvailable if available, otherwise false
      template <class Obj, bool defaultValue = false >
      struct SelectPreconditioning
      {
      private:
        template <typename T, typename = bool>
        struct CheckMember : public std::false_type { };

        template <typename T>
        struct CheckMember<T, decltype((void) T::preconditioningAvailable, true)> : public std::true_type { };

        template <class T, bool>
        struct SelectValue
        {
          static const bool value = defaultValue;
          typedef BaseType type;
        };

        template <class T>
        struct SelectValue< T, true >
        {
          static const bool value = T::preconditioningAvailable;
          typedef typename T::PreconditionerType type;
        };
      public:
        static constexpr bool value = SelectValue< Obj, CheckMember< Obj >::value >::value;
        typedef typename SelectValue< Obj, CheckMember< Obj >::value >::type type;
      };

    public:
      //! type of operator's Jacobian
      typedef JacobianOperator JacobianOperatorType;

      //! type of operator to invert
      typedef DifferentiableOperator< JacobianOperatorType > OperatorType;

      //! type of linear inverse operator
      typedef LInvOp LinearInverseOperatorType;

      //! type of preconditioner for linear solver
      static constexpr bool preconditioningAvailable = SelectPreconditioning< LinearInverseOperatorType > :: value;
      typedef typename SelectPreconditioning< LinearInverseOperatorType > :: type PreconditionerType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef NewtonParameter<typename LinearInverseOperatorType::SolverParameterType> ParameterType;

      typedef std::function< bool ( const RangeFunctionType &w, const RangeFunctionType &dw, double residualNorm ) > ErrorMeasureType;

      /** constructor
       *
       *  \param[in]  jInv       linear inverse operator (will be move constructed)
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */

      /** constructor
       *
       *  \param[in]  jInv        linear inverse operator (will be move constructed)
       *  \param[in]  epsilon     tolerance for norm of residual
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */

      // main constructor
      NewtonInverseOperator ( LinearInverseOperatorType jInv, const DomainFieldType &epsilon, const ParameterType &parameter )
        : verbose_( parameter.verbose() ),
          maxLineSearchIterations_( parameter.maxLineSearchIterations() ),
          jInv_( std::move( jInv ) ),
          parameter_(parameter),
          lsMethod_( parameter.lineSearch() ),
          finished_( [ epsilon ] ( const RangeFunctionType &w, const RangeFunctionType &dw, double res ) { return res < epsilon; } ),
          linearToleranceStrategy_ ( parameter.linearToleranceStrategy() ),
          eisenstatWalker_ ( epsilon )
      {
        if (linearToleranceStrategy_ == ParameterType::LinearToleranceStrategy::eisenstatwalker) {
          if (parameter_.linear().errorMeasure() != LinearSolver::ToleranceCriteria::residualReduction) {
            DUNE_THROW( InvalidStateException, "Parameter `newton.linear.errormeasure` selecting the tolerance criteria in the linear solver must be `residualreduction` when using Eisenstat-Walker." );
          }
        }
      }


      /** constructor
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */
      /*
      explicit NewtonInverseOperator ( const ParameterType &parameter )
        : NewtonInverseOperator( parameter.tolerance(), parameter )
      {}
      */

      explicit NewtonInverseOperator ( const ParameterType &parameter = ParameterType( Parameter::container() ) )
        : NewtonInverseOperator( parameter.tolerance(), parameter )
      {
        // std::cout << "in Newton inv op should use:" << parameter.linear().solverMethod({SolverParameter::gmres,SolverParameter::bicgstab,SolverParameter::minres}) << std::endl;
      }

      /** constructor
       *
       *  \param[in]  epsilon  tolerance for norm of residual
       */
      NewtonInverseOperator ( const DomainFieldType &epsilon, const ParameterType &parameter )
        : NewtonInverseOperator(
            LinearInverseOperatorType( parameter.linear() ),
            epsilon, parameter )
      {}

      NewtonInverseOperator ( const DomainFieldType &epsilon,
                              const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( epsilon, ParameterType( parameter ) )
      {}


      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */

      void setErrorMeasure ( ErrorMeasureType finished ) { finished_ = std::move( finished ); }

      EisenstatWalkerStrategy& eisenstatWalker () { return eisenstatWalker_; }

      void bind ( const OperatorType &op ) { op_ = &op; }

      void bind ( const OperatorType &op, const PreconditionerType& preconditioner )
      {
        bind( op );
        if( preconditioningAvailable )
          preconditioner_ = &preconditioner;
        else
          std::cerr << "WARNING: linear inverse operator does not support external preconditioning!" << std::endl;
      }

      void unbind () { op_ = nullptr; preconditioner_ = nullptr; }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const { return iterations_; }
      void setMaxIterations ( int maxIterations ) { parameter_.setMaxIterations( maxIterations ); }
      int linearIterations () const { return linearIterations_; }
      void setMaxLinearIterations ( int maxLinearIterations ) { parameter_.setMaxLinearIterations( maxLinearIterations ); }
      void updateLinearTolerance () const {
        if (linearToleranceStrategy_ == ParameterType::LinearToleranceStrategy::eisenstatwalker) {
          double newTol = eisenstatWalker_.nextLinearTolerance( delta_ );
          jInv_.parameter().setTolerance( newTol );
        }
      }
      bool verbose() const { return verbose_ && Parameter::verbose( Parameter::solverStatistics ); }
      double residual () const { return delta_; }

      NewtonFailure failed () const
      {
        // check for finite |residual| - this also works for -ffinite-math-only (gcc)
        // nan test not working with optimization flags...
        if( !(delta_ < std::numeric_limits< DomainFieldType >::max()) || std::isnan( delta_ ) )
          return NewtonFailure::InvalidResidual;
        else if( iterations_ >= parameter_.maxIterations() )
          return NewtonFailure::TooManyIterations;
        else if( linearIterations_ >= parameter_.maxLinearIterations() )
          return NewtonFailure::TooManyLinearIterations;
        else if( linearIterations_ < 0)
          return NewtonFailure::LinearSolverFailed;
        else if( !stepCompleted_ )
          return NewtonFailure::LineSearchFailed;
        else
          return NewtonFailure::Success;
      }

      bool converged () const { return failed() == NewtonFailure::Success; }

      virtual int lineSearch(RangeFunctionType &w, RangeFunctionType &dw,
                   const DomainFunctionType &u, DomainFunctionType &residual) const
      {
        double deltaOld = delta_;
        delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
        if (lsMethod_ == ParameterType::LineSearchMethod::none)
          return 0;
        if (failed() == NewtonFailure::InvalidResidual)
        {
          double test = dw.scalarProductDofs( dw );
          if (! (test < std::numeric_limits< DomainFieldType >::max() &&
                 !std::isnan(test)) )
            delta_ = 2.*deltaOld; // enter line search
        }
        double factor = 1.0;
        int noLineSearch = (delta_ < deltaOld)?1:0;
        int lineSearchIteration = 0;
        const bool lsVerbose = verbose() && Parameter::verbose( Parameter::extendedStatistics );
        while (delta_ >= deltaOld)
        {
          double deltaPrev = delta_;
          factor *= 0.5;
          if( lsVerbose )
            std::cout << "    line search:" << delta_ << ">" << deltaOld << std::endl;
          if (std::abs(delta_-deltaOld) < 1e-5*delta_) // || !converged()) // line search not working
            return -1;  // failed
          w.axpy(factor,dw);
          (*op_)( w, residual );
          residual -= u;
          delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
          if (std::abs(delta_-deltaPrev) < 1e-15)
            return -1;
          if (failed() == NewtonFailure::InvalidResidual)
            delta_ = 2.*deltaOld; // remain in line search

          ++lineSearchIteration;
          if( lineSearchIteration >= maxLineSearchIterations_ )
            return -1; // failed
        }
        return noLineSearch;
      }

    protected:
      void bindOperatorAndPreconditioner( JacobianOperatorType& jOp ) const
      {
        // if preconditioner was given pass it on to linear solver
        if constexpr ( preconditioningAvailable )
        {
          if( preconditioner_ )
          {
            jInv_.bind( jOp, *preconditioner_ );
            return ;
          }
        }

        // if preconditioning not enabled or set, then only set jOp
        jInv_.bind( jOp );
      }

      // hold pointer to jacobian operator, if memory reallocation is needed, the operator should know how to handle this.
      template< class ... Args>
      JacobianOperatorType& jacobian ( Args && ... args ) const
      {
        if( !jOp_ )
          jOp_.reset( new JacobianOperatorType( std::forward< Args >( args ) ...) ); //, parameter_.parameter() ) );
        return *jOp_;
      }

    protected:
      const OperatorType *op_ = nullptr;
      const PreconditionerType* preconditioner_ = nullptr;

      const bool verbose_;
      const int maxLineSearchIterations_;

      mutable DomainFieldType delta_;
      mutable int iterations_;
      mutable int linearIterations_;
      mutable LinearInverseOperatorType jInv_;
      mutable std::unique_ptr< JacobianOperatorType > jOp_;
      ParameterType parameter_;
      mutable int stepCompleted_;
      typename ParameterType::LineSearchMethod lsMethod_;
      ErrorMeasureType finished_;
      typename ParameterType::LinearToleranceStrategy linearToleranceStrategy_;
      EisenstatWalkerStrategy eisenstatWalker_;
    };


    // Implementation of NewtonInverseOperator
    // ---------------------------------------

    template< class JacobianOperator, class LInvOp >
    inline void NewtonInverseOperator< JacobianOperator, LInvOp >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      assert( op_ );

      DomainFunctionType residual( u );
      RangeFunctionType dw( w );
      JacobianOperatorType& jOp = jacobian( "jacobianOperator", dw.space(), u.space(), parameter_.solverParameter() );

      stepCompleted_ = true;
      iterations_ = 0;
      linearIterations_ = 0;
      // compute initial residual
      (*op_)( w, residual );
      residual -= u;
      delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
      updateLinearTolerance();

      const bool newtonVerbose = verbose();
      if( newtonVerbose )
        std::cout << "Newton iteration " << iterations_ << ": |residual| = " << delta_;
      while( true )
      {
        if( newtonVerbose )
          std::cout << std::endl;
        // evaluate operator's Jacobian
        (*op_).jacobian( w, jOp );

        // David: With this factor, the tolerance of CGInverseOp is the absolute
        //        rather than the relative error
        //        (see also dune-fem/dune/fem/solver/krylovinverseoperators.hh)

        // bind jOp to jInv including preconditioner if enabled and set
        bindOperatorAndPreconditioner( jOp );

        if ( parameter_.maxLinearIterations() - linearIterations_ <= 0 )
          break;
        jInv_.setMaxIterations( parameter_.maxLinearIterations() - linearIterations_ );

        dw.clear();
        jInv_( residual, dw );
        if (jInv_.iterations() < 0) // iterations are negative if solver didn't converge
        {
          linearIterations_ = jInv_.iterations();
          break;
        }
        linearIterations_ += jInv_.iterations();
        w -= dw;

        (*op_)( w, residual );
        residual -= u;
        int ls = lineSearch(w,dw,u,residual);
        stepCompleted_ = ls >= 0;
        updateLinearTolerance();
        ++iterations_;
        if( newtonVerbose )
          std::cout << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::flush;
        // if ( (ls==1 && finished_(w, dw, delta_)) || !converged())
        if ( (finished_(w, dw, delta_)) || !converged())
        {
          if( newtonVerbose )
          {
            std::cout << std::endl;
            std::cout << "Linear iterations: " << linearIterations_ << std::endl;
          }
          break;
        }
      }
      if( newtonVerbose )
        std::cout << std::endl;

      jInv_.unbind();
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
