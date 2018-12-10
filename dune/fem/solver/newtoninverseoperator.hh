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

    // NewtonParameter
    // ---------------

    struct NewtonParameter
      : public Dune::Fem::LocalParameter< NewtonParameter, NewtonParameter >
    {
      protected:

      std::shared_ptr<SolverParameter> baseParam_;
      // key prefix, default is fem.solver.newton. (can be overloaded by user)
      const std::string keyPrefix_;

      ParameterReader parameter_;

    public:
      NewtonParameter( const SolverParameter& baseParameter, const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( baseParameter.clone() ),
          keyPrefix_( keyPrefix ),
          parameter_( parameter )
      {}

      explicit NewtonParameter( const SolverParameter& baseParameter, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( baseParameter.clone() ),
          keyPrefix_( "fem.solver.newton." ),
          parameter_( parameter )
      {}

      NewtonParameter( const ParameterReader &parameter = Parameter::container() )
        : baseParam_( std::make_shared<SolverParameter>("fem.solver.newton.linear.", parameter) ),
          keyPrefix_( "fem.solver.newton." ),
          parameter_( parameter )
      {}

      NewtonParameter( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( std::make_shared<SolverParameter>(keyPrefix + "linear.", parameter) ),
          keyPrefix_( keyPrefix ),
          parameter_( parameter )
      {}

      const ParameterReader &parameter () const { return parameter_; }
      const SolverParameter& solverParameter () const { return *baseParam_; }

      const SolverParameter& linear () const { return *baseParam_; }

      //These methods affect the nonlinear solver
      virtual double tolerance () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "tolerance", 1e-6 );
      }

      virtual void setTolerance ( const double tol )
      {
        Parameter::append( keyPrefix_ + "tolerance", std::to_string(tol), true );
      }

      virtual bool verbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >(keyPrefix_ +  "verbose", v );
      }

      virtual int maxIterations () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxiterations", std::numeric_limits< int >::max() );
      }

      virtual void setMaxIterations ( const int maxIter )
      {
        Parameter::append( keyPrefix_ + "maxiterations", std::to_string(maxIter), true);
      }

      //Maximum Linear Iterations in total
      //!= max iterations of each linear solve
      virtual int maxLinearIterations () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
      }

      virtual void setMaxLinearIterations ( const int maxLinearIter )
      {
        Parameter::append( keyPrefix_ + "maxlineariterations", std::to_string(maxLinearIter), true);
      }

      virtual int maxLineSearchIterations () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlinesearchiterations", std::numeric_limits< int >::max() );
      }

      virtual void setMaxLineSearchIterations ( const int maxLineSearchIter )
      {
        Parameter::append( keyPrefix_ + "maxlinesearchiterations", std::to_string(maxLineSearchIter), true);
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

      //deprecated methods
      // split into removal of Parameter at the end
      // and forward to the linear parameters
      [[deprecated("Replaced by tolerance()")]]
      virtual double toleranceParameter () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "tolerance", 1e-6 );
      }

      [[deprecated("please use the linear solver parameters instead")]]
      virtual double linAbsTolParameter ( const double &tolerance )  const
      {
        if(parameter_.exists(keyPrefix_ + "linabstol"))
          return parameter_.getValue< double >(keyPrefix_ +  "linabstol", tolerance / 8 );
        return linear().linAbsTol();
      }

      [[deprecated("please use the linear solver parameters instead")]]
      virtual double linReductionParameter ( const double &tolerance ) const
      {
        if(parameter_.exists(keyPrefix_ + "linreduction"))
          return parameter_.getValue< double >( keyPrefix_ + "linreduction", tolerance / 8 );
        return linear().linReduction();
      }

      [[deprecated("Replaced by verbose ()")]]
      virtual bool newtonVerbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >(keyPrefix_ +  "verbose", v );
      }

      [[deprecated("please use the linear solver parameters instead")]]
      virtual bool linearSolverVerbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >( keyPrefix_ + "linear.verbose", v );
      }

      [[deprecated("Replaced by maxIterations ()")]]
      virtual int maxIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxiterations", std::numeric_limits< int >::max() );
      }

      [[deprecated("Replaced by maxLinearIterations")]]
      virtual int maxLinearIterationsParameter () const
      {
        if(parameter_.exists(keyPrefix_ + "maxlineariterations"))
          return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
        return linear().maxLinearIterations();
      }

      [[deprecated]]
      virtual int maxLineSearchIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlinesearchiterations", std::numeric_limits< int >::max() );
      }

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
      TooManyLinearIterations = 6
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

    public:
      //! type of operator's Jacobian
      typedef JacobianOperator JacobianOperatorType;

      //! type of operator to invert
      typedef DifferentiableOperator< JacobianOperatorType > OperatorType;

      //! type of linear inverse operator
      typedef LInvOp LinearInverseOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      typedef NewtonParameter ParametersType;

      typedef std::function< bool ( const RangeFunctionType &w, const RangeFunctionType &dw, double residualNorm ) > ErrorMeasureType;

      /** constructor
       *
       *  \param[in]  jInv       linear inverse operator (will be move constructed)
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */

      NewtonInverseOperator ( LinearInverseOperatorType jInv, const NewtonParameter &parameter )
        : NewtonInverseOperator( std::move( jInv ), parameter.tolerance(), parameter )
      {}

      explicit NewtonInverseOperator ( LinearInverseOperatorType jInv,
                                       const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( std::move( jInv ), ParametersType( parameter ) )
      {}

      /** constructor
       *
       *  \param[in]  jInv        linear inverse operator (will be move constructed)
       *  \param[in]  epsilon     tolerance for norm of residual
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */

      NewtonInverseOperator ( LinearInverseOperatorType jInv, const DomainFieldType &epsilon, const NewtonParameter &parameter )
        : verbose_( parameter.verbose() && MPIManager::rank () == 0 ),
          maxLineSearchIterations_( parameter.maxLineSearchIterations() ),
          jInv_( std::move( jInv ) ),
          parameter_(parameter),
          lsMethod_( parameter.lineSearch() ),
          finished_( [ epsilon ] ( const RangeFunctionType &w, const RangeFunctionType &dw, double res ) { return res < epsilon; } )
      {}

      NewtonInverseOperator ( LinearInverseOperatorType jInv, const DomainFieldType &epsilon,
                              const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( std::move( jInv ), epsilon, ParametersType( parameter ) )
      {}


      /** constructor
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */
      explicit NewtonInverseOperator ( const NewtonParameter &parameter )
        : NewtonInverseOperator( parameter.tolerance(), parameter )
      {}

      explicit NewtonInverseOperator ( const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( ParametersType( parameter ) )
      {}

      /** constructor
       *
       *  \param[in]  epsilon  tolerance for norm of residual
       */
      NewtonInverseOperator ( const DomainFieldType &epsilon, const NewtonParameter &parameter )
        : NewtonInverseOperator(
            LinearInverseOperatorType( parameter.linear().linReduction( epsilon ),
                                       parameter.linear().linAbsTol( epsilon ),
                                       parameter.linear().maxLinearIterations(),
                                       parameter.linear().verbose(),
                                       parameter.linear() ),
            epsilon, parameter )
      {}

      NewtonInverseOperator ( const DomainFieldType &epsilon,
                              const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( epsilon, ParametersType( parameter ) )
      {}


      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */


      NewtonInverseOperator ( const OperatorType &op, const NewtonParameter &parameter )
        : NewtonInverseOperator( parameter )
      {
        bind( op );
      }

      explicit NewtonInverseOperator ( const OperatorType &op,
                                       const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( parameter )
      {
        bind( op );
      }

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  epsilon  tolerance for norm of residual
       */

      NewtonInverseOperator ( const OperatorType &op, const DomainFieldType &epsilon,
                              const NewtonParameter &parameter )
        : NewtonInverseOperator( epsilon, parameter )
      {
        bind( op );
      }

      NewtonInverseOperator ( const OperatorType &op, const DomainFieldType &epsilon,
                              const ParameterReader &parameter = Parameter::container() )
        : NewtonInverseOperator( epsilon, parameter )
      {
        bind( op );
      }

      void setErrorMeasure ( ErrorMeasureType finished ) { finished_ = std::move( finished ); }

      void bind ( const OperatorType &op ) { op_ = &op; }

      void unbind () { op_ = nullptr; }

      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const { return iterations_; }
      void setMaxIterations ( int maxIterations ) { parameter_.setMaxIterations( maxIterations ); }
      int linearIterations () const { return linearIterations_; }
      void setMaxLinearIterations ( int maxLinearIterations ) { parameter_.setMaxLinearIterations( maxLinearIterations ); }
      bool verbose() const { return verbose_; }

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
        if (lsMethod_ == NewtonParameter::LineSearchMethod::none)
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
        while (delta_ >= deltaOld)
        {
          double deltaPrev = delta_;
          factor *= 0.5;
          if( verbose() )
            std::cerr << "    line search:" << delta_ << ">" << deltaOld << std::endl;
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
      // hold pointer to jacobian operator, if memory reallocation is needed, the operator should know how to handle this.
      template< class ... Args>
      JacobianOperatorType& jacobian ( Args && ... args ) const
      {
        if( !jOp_ )
          jOp_.reset( new JacobianOperatorType( std::forward< Args >( args ) ...) ); //, parameter_.parameter() ) );
        return *jOp_;
      }

    private:
      const OperatorType *op_ = nullptr;

      const bool verbose_;
      const int maxLineSearchIterations_;

      mutable DomainFieldType delta_;
      mutable int iterations_;
      mutable int linearIterations_;
      mutable LinearInverseOperatorType jInv_;
      mutable std::unique_ptr< JacobianOperatorType > jOp_;
      NewtonParameter parameter_;
      mutable int stepCompleted_;
      NewtonParameter::LineSearchMethod lsMethod_;
      ErrorMeasureType finished_;
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
      JacobianOperatorType& jOp = jacobian( "jacobianOperator", dw.space(), u.space()); // , parameter_.parameter() );

      stepCompleted_ = true;
      iterations_ = 0;
      linearIterations_ = 0;
      // compute initial residual
      (*op_)( w, residual );
      residual -= u;
      delta_ = std::sqrt( residual.scalarProductDofs( residual ) );

      if( verbose() )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_;
      while( true )
      {
        if( verbose() )
          std::cerr << std::endl;
        // evaluate operator's jacobian
        (*op_).jacobian( w, jOp );

        // David: With this factor, the tolerance of CGInverseOp is the absolute
        //        rather than the relative error
        //        (see also dune-fem/dune/fem/solver/krylovinverseoperators.hh)
        jInv_.bind( jOp );
        jInv_.setMaxIterations( parameter_.maxLinearIterations() - linearIterations_ );

        dw.clear();
        jInv_( residual, dw );
        linearIterations_ += jInv_.iterations();
        w -= dw;

        (*op_)( w, residual );
        residual -= u;
        int ls = lineSearch(w,dw,u,residual);
        stepCompleted_ = ls >= 0;
        ++iterations_;
        if( verbose() )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::flush;
        // if ( (ls==1 && finished_(w, dw, delta_)) || !converged())
        if ( (finished_(w, dw, delta_)) || !converged())
        {
          if( verbose() )
            std::cerr << std::endl;
          break;
        }
      }
      if( verbose() )
        std::cerr << std::endl;

      jInv_.unbind();
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
