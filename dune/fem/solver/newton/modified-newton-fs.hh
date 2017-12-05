#ifndef DUNE_FEM_MODIFIED_NEWTON_FS_HH
#define DUNE_FEM_MODIFIED_NEWTON_FS_HH

#include <dune/fem/solver/newtoninverseoperator.hh>

namespace Dune {

  namespace Fem {

    /**
    *  solve nonlinear systems by a Newton method with step size control
    *  Reimplementation of newton_fs.c in alberta_util/src of ALBERTA
    */
    template<class JacobianOperator, class LInvOp>
    class FSNewtonInverseOperator
      : public Dune::Fem::NewtonInverseOperator<JacobianOperator, LInvOp>
    {
      typedef Dune::Fem::NewtonInverseOperator<JacobianOperator, LInvOp> BaseType;
    public:
      using BaseType::NewtonInverseOperator;
      typedef JacobianOperator JacobianOperatorType;
      using BaseType::converged;
      using BaseType::jacobian;
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;
      typedef typename BaseType::LinearInverseOperatorType LinearInverseOperatorType;

      using BaseType::op_;
      using BaseType::jInv_;

      using BaseType::tolerance_;
      using BaseType:: verbose_;
      using BaseType:: maxLinearIterations_;

      using BaseType::delta_;
      using BaseType::iterations_;
      using BaseType::linearIterations_;

      //max Number of controlling steps
      int mmax_ = Dune::Fem::Parameter::getValue<int>("fem.solver.newton.stepsizeiterations",10);

      void operator()(const DomainFunctionType& u, RangeFunctionType& w) const
      {
        //double stepSizeTol = 0.1 ;

        DomainFunctionType y(w);
        y.clear();
        double deltaY=1e308;
        double tau=1.;
        bool halved=true;


        DomainFunctionType residual;
        residual.clear();
        DomainFunctionType residualOfY;
        residualOfY.clear();
        RangeFunctionType dw;
        dw.clear();
        JacobianOperatorType& jOp = jacobian( "jacobianOperator", dw.space(), u.space() );


        // compute initial residual
        (*op_)( w, residual );
        residual -= u;

        delta_ = std::sqrt( residual.scalarProductDofs( residual ) );

        for( iterations_ = 0, linearIterations_ = 0; converged() && (delta_ > tolerance_); ++iterations_ )
        {

          if( verbose_ )
            std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;

          (*op_).jacobian( w, jOp );
          jInv_.bind( jOp );
          jInv_.setMaxIterations( maxLinearIterations_ - linearIterations_ );
          dw.clear();
          jInv_( residual, dw );
          linearIterations_ += jInv_.iterations();

          if (!halved)
          {
            tau = tau < 0.5 ? 2.0*tau : 1.0;
          }

          for (int j = 0; j <= mmax_; j++)
          {
            y.assign(w);
            y.axpy(-tau,dw);

            (*op_)(y, residualOfY);
            residualOfY -= u;
            deltaY = std::sqrt( residualOfY.scalarProductDofs( residualOfY ) );

            /*--- aim: |F(u_k+\tau d)| \le (1-0.5\tau) |F(u)| --------------------------*/
            if (deltaY <= (1.0 - 0.5*tau)*delta_)
            {
              halved = false;
              break;
            }
            else
            {
              halved=true;
              tau *= 0.5;
            }
          }
          w.assign(y);
          residual.assign(residualOfY);
          delta_=deltaY;
        }
        jInv_.unbind();
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;
      }   // operator ()

    };   //class MyNewtonInverseOperator
  }   // namespace Fem

}   // namespace Dune



#endif // DUNE_FEM_MODIFIED_NEWTON_FS_HH
