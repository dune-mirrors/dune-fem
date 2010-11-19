// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_NEWTONINVERSEOPERATOR_HH
#define DUNE_FEM_NEWTONINVERSEOPERATOR_HH

#include <cfloat>
#include <iostream>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>

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
    template< class Op, class LInvOp >
    class NewtonInverseOperator 
    : public Operator< typename Op::RangeFunctionType, typename Op::DomainFunctionType > 
    {
      typedef NewtonInverseOperator< Op, LInvOp > ThisType;
      typedef Operator< typename Op::RangeFunctionType, typename Op::DomainFunctionType > BaseType;

    public:
      //! type of operator to invert
      typedef Op OperatorType;

      //! type of operator's Jacobian
      typedef typename OperatorType::JacobianOperatorType JacobianOperatorType;

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
      const DomainFieldType tolerance_;
      const bool verbose_;
      const int maxIterations_;
      const int maxLinearIterations_;

      mutable int iterations_;
      mutable int linearIterations_;
    };


    
    template< class Op, class LInvOp >
    inline void NewtonInverseOperator< Op, LInvOp >
      ::operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      const DomainFieldType reduction = tolerance_ / 8;
      const DomainFieldType absLimit = tolerance_ / 8;

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
        const LinearInverseOperatorType jInv( jOp, reduction, absLimit / delta );
        
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
