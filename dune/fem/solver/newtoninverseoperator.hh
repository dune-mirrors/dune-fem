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
      // key prefix, default is fem.ode.newton. (can be overloaded by user)
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
        : baseParam_( std::make_shared<SolverParameter>(parameter) ),
          keyPrefix_( "fem.solver.newton." ),
          parameter_( parameter )
      {}

      NewtonParameter( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : baseParam_( std::make_shared<SolverParameter>(keyPrefix, parameter) ),
          keyPrefix_( keyPrefix ),
          parameter_( parameter )
      {}

      const ParameterReader &parameter () const { return parameter_; }
      const SolverParameter& solverParameter () const { return *baseParam_; }

      virtual double toleranceParameter () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "tolerance", 1e-6 );
      }

      virtual double linAbsTolParameter ( const double &tolerance )  const
      {
        return parameter_.getValue< double >(keyPrefix_ +  "linabstol", tolerance / 8 );
      }

      virtual double linReductionParameter ( const double &tolerance ) const
      {
        return parameter_.getValue< double >( keyPrefix_ + "linreduction", tolerance / 8 );
      }

      virtual bool newtonVerbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >(keyPrefix_ +  "verbose", v );
      }

      virtual bool linearSolverVerbose () const
      {
        const bool v = baseParam_? baseParam_->verbose() : false;
        return parameter_.getValue< bool >( keyPrefix_ + "linear.verbose", v );
      }

      virtual int maxIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxiterations", std::numeric_limits< int >::max() );
      }

      virtual int maxLinearIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
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
        : NewtonInverseOperator( std::move( jInv ), parameter.toleranceParameter(), parameter )
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
        : verbose_( parameter.newtonVerbose() && MPIManager::rank () == 0 ),
          maxIterations_( parameter.maxIterationsParameter() ),
          maxLinearIterations_( parameter.maxLinearIterationsParameter() ),
          maxLineSearchIterations_( parameter.maxLineSearchIterationsParameter() ),
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
        : NewtonInverseOperator( parameter.toleranceParameter(), parameter )
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
            LinearInverseOperatorType( parameter.linReductionParameter( epsilon ),
                                       parameter.linAbsTolParameter( epsilon ),
                                       parameter.maxLinearIterationsParameter(),
                                       parameter.linearSolverVerbose(),
                                       parameter.parameter() ),
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
      void setMaxIterations ( int maxIterations ) { maxIterations_ = maxIterations; }
      int linearIterations () const { return linearIterations_; }
      void setMaxLinearIterations ( int maxLinearIterations ) { maxLinearIterations_ = maxLinearIterations; }
      bool verbose() const { return verbose_; }

      NewtonFailure failed () const
      {
        // check for finite |residual| - this also works for -ffinite-math-only (gcc)
        // nan test not working with optimization flags...
        if( !(delta_ < std::numeric_limits< DomainFieldType >::max()) || std::isnan( delta_ ) )
          return NewtonFailure::InvalidResidual;
        else if( iterations_ >= maxIterations_ )
          return NewtonFailure::TooManyIterations;
        else if( linearIterations_ >= maxLinearIterations_ )
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
      const int maxIterations_;
      const int maxLinearIterations_;
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
        jInv_.setMaxIterations( maxLinearIterations_ - linearIterations_ );

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
