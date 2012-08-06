// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH

#include <cfloat>
#include <iostream>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/differentiableoperator.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class NewtonInverseOperator >
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
    : public Fem :: Operator< typename JacobianOperator :: DomainFunctionType, 
                              typename JacobianOperator :: RangeFunctionType > 
    {
      typedef NewtonInverseOperator< JacobianOperator, LInvOp > ThisType;
      typedef Fem::Operator< typename JacobianOperator :: DomainFunctionType,
                             typename JacobianOperator :: RangeFunctionType > BaseType;
    public:
      //! type of operator's Jacobian
      typedef JacobianOperator JacobianOperatorType;
      
      //! type of operator to invert
      typedef DifferentiableOperator<JacobianOperatorType> OperatorType;

      //! type of linear inverse operator
      typedef LInvOp LinearInverseOperatorType;

      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *
       *  \note The tolerance is read from the paramter
       *        <b>fem.solver.newton.tolerance</b>
       */
      explicit NewtonInverseOperator ( OperatorType &op )
      : op_( op ),
        tolerance_( toleranceParameter() ),
        linAbsTol_( linAbsTolParameter( tolerance_ ) ),
        linReduction_( linReductionParameter( tolerance_ ) ),
        verbose_( verbosityParameter() ),
        maxIterations_( maxIterationsParameter() ),
        maxLinearIterations_( maxLinearIterationsParameter() )
      {}

      /** constructor
       *
       *  \param[in]  op       operator to invert
       *  \param[in]  epsilon  tolerance for norm of residual
       */
      NewtonInverseOperator ( OperatorType &op, const DomainFieldType &epsilon )
      : op_( op ),
        tolerance_( epsilon ),
        linAbsTol_( linAbsTolParameter( tolerance_ ) ),
        linReduction_( linReductionParameter( tolerance_ ) ),
        verbose_( verbosityParameter() ),
        maxIterations_( maxIterationsParameter() ),
        maxLinearIterations_( maxLinearIterationsParameter() )
      {}
        
      virtual void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const;

      int iterations () const
      {
        return iterations_;
      }

      int linearIterations () const
      {
        return linearIterations_;
      }

      bool converged () const
      {
        return (iterations_ < maxIterations_) && (linearIterations_ < maxLinearIterations_);
      }

    private:
      static DomainFieldType toleranceParameter ()
      {
        return Dune::Parameter::getValue< DomainFieldType >( "fem.solver.newton.tolerance", 1e-6 );
      }

      static DomainFieldType linAbsTolParameter ( const DomainFieldType &tolerance )
      {
        return Dune::Parameter::getValue< DomainFieldType >( "fem.solver.newton.linabstol", tolerance / 8 );
      }

      static DomainFieldType linReductionParameter ( const DomainFieldType &tolerance )
      {
        return Dune::Parameter::getValue< DomainFieldType >( "fem.solver.newton.linreduction", tolerance / 8 );
      }

      static bool verbosityParameter ()
      {
        const bool v = Parameter::getValue< bool >( "fem.solver.verbose", false );
        return Parameter::getValue< bool >( "fem.solver.newton.verbose", v );
      }

      static int maxIterationsParameter ()
      {
        return Parameter::getValue< int >( "fem.solver.newton.maxiterations", std::numeric_limits< int >::max() );
      }

      static int maxLinearIterationsParameter ()
      {
        return Parameter::getValue< int >( "fem.solver.newton.maxlineariterations", std::numeric_limits< int >::max() );
      }

      OperatorType &op_;
      const DomainFieldType tolerance_, linAbsTol_, linReduction_;;
      const bool verbose_;
      const int maxIterations_;
      const int maxLinearIterations_;

      mutable int iterations_;
      mutable int linearIterations_;
    };


    
    template< class JacobianOperator, class LInvOp >
    inline void NewtonInverseOperator< JacobianOperator, LInvOp >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      DomainFunctionType residual( u );
      RangeFunctionType dw( w );
      JacobianOperatorType jOp( "jacobianOperator", dw.space(), u.space() );

      // compute initial residual
      op_( w, residual );
      residual -= u;
      DomainFieldType delta = std::sqrt( residual.scalarProductDofs( residual ) );

      for( iterations_ = 0, linearIterations_ = 0; converged() && (delta > tolerance_); ++iterations_ )
      {
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta << std::endl;

        // evaluate operator's jacobian
        op_.jacobian( w, jOp );
        
        // David: With this factor, the tolerance of CGInverseOp is the absolute
        //        rather than the relative error
        //        (see also dune-fem/dune/fem/solver/inverseoperators.hh)
        const int remLinearIts = maxLinearIterations_ - linearIterations_;
        const LinearInverseOperatorType jInv( jOp, linReduction_, linAbsTol_ / delta, remLinearIts );
        
        dw.clear();
        jInv( residual, dw );
        linearIterations_ += jInv.iterations();
        w -= dw;
        
        op_( w, residual );
        residual -= u;
        delta = std::sqrt( residual.scalarProductDofs( residual ) );
      }
      if( verbose_ )
        std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta << std::endl;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
