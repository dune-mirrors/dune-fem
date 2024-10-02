#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH

#include <cassert>
#include <cmath>

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <regex>
#include <utility>

#include <dune/common/timer.hh>
#include <dune/common/exceptions.hh>

#include <dune/fem/common/staticlistofint.hh>
#include <dune/fem/solver/parameter.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class EisenstatWalkerStrategy
     *
     * \brief Adaptive tolerance selection for linear solver.
     *
     * \note Prevents over-solving linear systems far away from solution.
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
      // key prefix, default is fem.solver.nonlinear. (can be overloaded by user)
      const std::string keyPrefix_;

      ParameterReader parameter_;

      void checkDeprecatedParameters() const
      {
        const std::string newton("newton.");
        const std::size_t pos = keyPrefix_.find( newton );
        if( pos != std::string::npos )
        {
          DUNE_THROW(InvalidStateException,"Keyprefix 'newton' is deprecated, replace with 'nonlinear'!");
        }

        const std::string params[]
          = { "tolerance", "lineSearch", "maxterations", "linear", "maxlinesearchiterations" };
        for( const auto& p : params )
        {
          std::string key( "fem.solver.newton." );
          key += p;
          if( parameter_.exists( key ) )
            DUNE_THROW(InvalidStateException,"Keyprefix 'newton' is deprecated, replace with 'nonlinear'!");
        }
      }

      std::string replaceNonLinearWithLinear( const std::string& keyPrefix ) const
      {
        if( keyPrefix.find( "nonlinear" ) != std::string::npos )
        {
          return std::regex_replace(keyPrefix, std::regex("nonlinear"), "linear" );
        }
        else
          return keyPrefix;
      }
    public:
      NewtonParameter( const SolverParam& baseParameter, const std::string keyPrefix = "fem.solver.nonlinear." )
        : baseParam_( static_cast< SolverParam* > (baseParameter.clone()) ),
          keyPrefix_( keyPrefix ),
          parameter_( baseParameter.parameter() )
      {
        checkDeprecatedParameters();
        checkForcingErrorMeasure();
      }

      template <class Parameter, std::enable_if_t<!std::is_base_of<SolverParam,Parameter>::value && !std::is_same<Parameter,ParameterReader>::value,int> i=0>
      NewtonParameter( const Parameter& solverParameter, const std::string keyPrefix = "fem.solver.nonlinear." )
        : baseParam_( new SolverParam(solverParameter) ),
          keyPrefix_( keyPrefix ),
          parameter_( solverParameter.parameter() )
      {
        checkDeprecatedParameters();
        checkForcingErrorMeasure();
      }

      template <class ParamReader, std::enable_if_t<!std::is_same<ParamReader,SolverParam>::value && std::is_same<ParamReader,ParameterReader>::value,int> i=0>
      NewtonParameter( const ParamReader &parameter, const std::string keyPrefix = "fem.solver.nonlinear." )
          // pass keyprefix for linear solvers, which is the same as keyprefix with nonlinear replaced by linear
        : baseParam_( std::make_shared<SolverParam>( replaceNonLinearWithLinear(keyPrefix), parameter) ),
          keyPrefix_( keyPrefix ),
          parameter_( parameter )
      {
        checkDeprecatedParameters();
        checkForcingErrorMeasure();
      }

      void checkForcingErrorMeasure()
      {
        if (forcing() == Forcing::eisenstatwalker)
        {
          baseParam_->setDefaultErrorMeasure(2);
          if (baseParam_->errorMeasure() != LinearSolver::ToleranceCriteria::residualReduction)
            DUNE_THROW( InvalidStateException, "Parameter `linear.errormeasure` selecting the tolerance criteria in the linear solver must be `residualreduction` when using Eisenstat-Walker." );
        }
      }

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

      // LineSearchMethod: none, simple
      LIST_OF_INT(LineSearchMethod,
                  none=0,
                  simple=1);

      virtual int lineSearch () const
      {
        if( parameter_.exists( keyPrefix_ + "lineSearch" ) )
        {
          std::cout << "WARNING: using old parameter name '" << keyPrefix_ + "lineSearch" << "',\n"
                    << "please switch to '" << keyPrefix_ + "linesearch" << "' (all lower caps)!" <<std::endl;
          return Forcing::to_id( parameter_.getEnum( keyPrefix_ + "lineSearch", LineSearchMethod::names(), LineSearchMethod::none ) );
        }
        return Forcing::to_id( parameter_.getEnum( keyPrefix_ + "linesearch", LineSearchMethod::names(), LineSearchMethod::none ) );
      }

      virtual void setLineSearch ( const int method )
      {
        Parameter::append( keyPrefix_ + "linesearch", LineSearchMethod::to_string(method), true );
      }

      // Forcing: none, eisenstatwalker
      LIST_OF_INT(Forcing,
                  none  =  0, // the provided linear solver tol is used in every iteration
                  eisenstatwalker=1); // Eistenstat-Walker criterion

      virtual int forcing () const
      {
        if( parameter_.exists( keyPrefix_ + "linear.tolerance.strategy" ) )
        {
          std::string keypref( keyPrefix_ );
          std::string femsolver("fem.solver.");
          size_t pos = keypref.find( femsolver );
          if (pos != std::string::npos)
          {
            // If found then erase it from string
            keypref.erase(pos, femsolver.length());
          }
          std::cout << "WARNING: using old parameter name '" << keypref + "linear.tolerance.strategy" << "',\n"
                    << "please switch to '" << keypref + "forcing" << "'!" <<std::endl;
          return Forcing::to_id( parameter_.getEnum( keyPrefix_ + "linear.tolerance.strategy", Forcing::names(), Forcing::none ) );
        }
        return Forcing::to_id( parameter_.getEnum( keyPrefix_ + "forcing", Forcing::names(), Forcing::none ) );
      }

      virtual void setForcing ( const int strategy )
      {
        Parameter::append( keyPrefix_ + "forcing", Forcing::to_string( strategy ), true );
      }

      //! return true if simplified Newton is to be used
      virtual bool simplified () const
      {
        return parameter_.getValue< bool >( keyPrefix_ + "simplified", 0 );
      }

      // allow to override the automatic choice of nonlinear or linear solver to
      // force nonlinear all the time
      virtual bool forceNonLinear () const
      {
        bool v = false;
        v = parameter_.getValue< bool >(keyPrefix_ +  "forcenonlinear", v );
        return v;
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
     *        paramter <b>fem.solver.nonlinear.verbose</b>; it defaults to
     *        <b>fem.solver.verbose</b>.
     *
     *  \note Similar to CG solver the initial guess should take the
     *        constraints into account. In this case we need the initial
     *        residual to be zero on the boundary, i.e., when calling
     *        operator()(u,w) then w=g+u should hold on the Dirichlet boundary.
     *        This way we do not explicitly need to call the constraints here.
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

      //! performance info about last solver call
      typedef Impl::SolverInfo SolverInfoType;

      /** constructor
       *
       *  \param[in]  jInv       linear inverse operator (will be move constructed)
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.nonlinear.tolerance</b>
       */

      /** constructor
       *
       *  \param[in]  jInv        linear inverse operator (will be move constructed)
       *  \param[in]  epsilon     tolerance for norm of residual
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.nonlinear.tolerance</b>
       */

      // main constructor
      NewtonInverseOperator ( LinearInverseOperatorType jInv, const DomainFieldType &epsilon, const ParameterType &parameter )
        : verbose_( parameter.verbose() ),
          maxLineSearchIterations_( parameter.maxLineSearchIterations() ),
          jInv_( std::move( jInv ) ),
          parameter_(parameter),
          lsMethod_( parameter.lineSearch() ),
          finished_( [ epsilon ] ( const RangeFunctionType &w, const RangeFunctionType &dw, double res ) { return res < epsilon; } ),
          forcing_ ( parameter.forcing() ),
          eisenstatWalker_ ( epsilon ),
          timing_(3, 0.0)
      {
      }


      /** constructor
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.nonlinear.tolerance</b>
       */
      explicit NewtonInverseOperator ( const ParameterType &parameter = ParameterType( Parameter::container() ) )
        : NewtonInverseOperator( parameter.tolerance(), parameter )
      {
        // std::cout << "in Newton inv op should use:" << parameter.linear().solverMethod({SolverParameter::gmres,SolverParameter::bicgstab,SolverParameter::minres}) << std::endl;
      }

      /** constructor
       *
       *  \param[in]  epsilon     tolerance for norm of residual
       *  \param[in]  parameter   parameter set for solver config.
       */
      NewtonInverseOperator ( const DomainFieldType &epsilon, const ParameterType &parameter )
        : NewtonInverseOperator(
            LinearInverseOperatorType( parameter.linear() ),
            epsilon, parameter )
      {}

      /** constructor
       *
       *  \param[in]  epsilon     tolerance for norm of residual
       *  \param[in]  parameter   parameter set for solver config.
       */
      NewtonInverseOperator ( const DomainFieldType &epsilon,
                              const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( epsilon, ParameterType( parameter ) )
      {}


      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.nonlinear.tolerance</b>
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
        if (forcing_ == ParameterType::Forcing::eisenstatwalker) {
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

      //! returns [overall, jacobian, solve] timings in seconds for last operator () call.
      const std::vector<double>& timing() const { return timing_; }

      //! return performance information about last solver run */
      SolverInfoType info() const { return SolverInfoType( converged(), linearIterations(), iterations(), timing() ); }

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
      const int lsMethod_;
      ErrorMeasureType finished_;
      const int forcing_;
      EisenstatWalkerStrategy eisenstatWalker_;

      mutable std::vector<double> timing_;
    };


    // Implementation of NewtonInverseOperator
    // ---------------------------------------

    template< class JacobianOperator, class LInvOp >
    inline void NewtonInverseOperator< JacobianOperator, LInvOp >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      assert( op_ );
      std::fill(timing_.begin(), timing_.end(), 0.0 );

      // obtain information about operator to invert
      const bool nonlinear = op_->nonlinear() || parameter_.forceNonLinear();

      Dune::Timer allTimer;
      DomainFunctionType residual( u );
      RangeFunctionType dw( w );

      Dune::Timer jacTimer;
      double jacTime = 0.0;
      JacobianOperatorType& jOp = jacobian( "jacobianOperator", dw.space(), u.space(), parameter_.solverParameter() );
      jacTime += jacTimer.elapsed();

      stepCompleted_ = true;
      iterations_ = 0;
      linearIterations_ = 0;
      // compute initial residual
      (*op_)( w, residual );   // r=S[w],   r=w-g on bnd
      residual -= u;           // r=S[w]-u, r=w-g-u on bnd (note: we should start with w=g+u on bnd so r=0)
      delta_ = std::sqrt( residual.scalarProductDofs( residual ) );
      updateLinearTolerance();

      // this is true for Newton, and false simplified Newton after first iteration
      bool computeJacobian = true;
      const bool simplifiedNewton = parameter_.simplified();

      const bool newtonVerbose = verbose() && nonlinear;
      if( newtonVerbose )
      {
        std::cout << "Start Newton: tol = " << parameter_.tolerance() << " (forcing = " << ParameterType::Forcing::to_string(forcing_) << " | linear tol = " << parameter_.linear().tolerance() << ")"<<std::endl;
        std::cout << "Newton iteration " << iterations_ << ": |residual| = " << delta_;
      }
      while( true )
      {
        if( newtonVerbose )
          std::cout << std::endl;

        if( computeJacobian )
        {
          // evaluate operator's Jacobian
          jacTimer.reset();
          (*op_).jacobian( w, jOp );

          // if simplified Newton is activated do not compute Jacobian again
          computeJacobian = ! simplifiedNewton;
        }

        // bind jOp to jInv including preconditioner if enabled and set
        bindOperatorAndPreconditioner( jOp );
        jacTime += jacTimer.elapsed();

        if ( parameter_.maxLinearIterations() - linearIterations_ <= 0 )
          break;
        jInv_.setMaxIterations( parameter_.maxLinearIterations() - linearIterations_ );

        dw.clear();
        jInv_( residual, dw );  // dw = DS[w]^{-1}(S[w]-u)
                                // dw = w-g-u on bnd
        if (jInv_.iterations() < 0) // iterations are negative if solver didn't converge
        {
          linearIterations_ = jInv_.iterations();
          break;
        }
        linearIterations_ += jInv_.iterations();
        w -= dw;                // w = w - DS[w]^{-1}(S[w]-u)
                                // w = g+u

        // the following only for nonlinear problems
        if( nonlinear )
        {
          // compute new residual
          (*op_)( w, residual ); // res = S[w]
                                 // res = w-g = g+u-g = u
          residual -= u;         // res = S[w] - u
                                 // res = 0

          // line search if enabled
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
        else // in linear case do not continue
          break ;
      } // end Newton loop

      if( newtonVerbose )
        std::cout << std::endl;

      jInv_.unbind();

      // store time measurements
      timing_[0] = allTimer.elapsed();
      timing_[1] = jacTime;
      timing_[2] = timing_[0] - jacTime;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
