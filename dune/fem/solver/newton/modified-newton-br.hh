#ifndef DUNE_FEM_MODIFIED_NEWTON_HH
#define DUNE_FEM_MODIFIED_NEWTON_HH

#include <dune/fem/solver/newtoninverseoperator.hh>

namespace Dune {

  namespace Fem {

    /**
    *  solve nonlinear systems by a Newton method with step size control
    *  by Bank and Rose, Numer. Math. 37 (1981) pp. 279-295.
    *  Reimplementation of newton_br.c in alberta_util/src of ALBERTA
    */

    template<class JacobianOperator, class LInvOp>
    class BRNewtonInverseOperator
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
      //stepsize control tolerance
      double stepSizeTol_ = Dune::Fem::Parameter::getValue<double>("fem.solver.newton.stepsizetol",0.1);

      void operator()(const DomainFunctionType& u, RangeFunctionType& w) const
      {
        DomainFunctionType y(w);
        y.clear();
        //residual of helper function
        double deltaY = 1e308;
        //stepsize
        double tau = 1.;
        // constant to control stepsize - adjusted in the loop
        double K = 0.0;


        DomainFunctionType residual;
        residual.clear();
        DomainFunctionType residualOfY;
        residualOfY.clear();
        //descent direction
        RangeFunctionType dw;
        dw.clear();
        JacobianOperatorType& jOp = jacobian( "jacobianOperator", dw.space(), u.space() );

        // compute initial residual
        (*op_)( w, residual );
        residual -= u;
        delta_=std::sqrt( residual.scalarProductDofs( residual ) );
        for( iterations_ = 0, linearIterations_ = 0; converged() && (delta_ > tolerance_); ++iterations_ )
        {
          if( verbose_ )
            std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;

          // update DF(x)
          (*op_).jacobian( w, jOp );
          jInv_.bind( jOp );
          jInv_.setMaxIterations( maxLinearIterations_ - linearIterations_ );
          dw.clear();
          jInv_( residual, dw );
          linearIterations_ += jInv_.iterations();

          for (int m = 0; m <= mmax_; m++)
          {
            /*--- aim: |F(u_k+\tau d)| \le ????--------------------------*/
            tau = 1.0 / (1.0 + K*delta_);
            //set y = w - tau*dw
            y.assign(w);
            y.axpy(-tau,dw);

            //calculate residual F(y) -u
            (*op_)(y, residualOfY);
            residualOfY -= u;
            deltaY = std::sqrt( residualOfY.scalarProductDofs( residualOfY ) );

            if ((1.0 - deltaY/delta_)/tau < stepSizeTol_)
            {
              K = (K == 0.0) ? 1.0 : 10.0*K;
            }
            else
            {
              K *= 0.1;
              break;
            }
          }//Stepsize Loop

          w.assign(y);
          residual.assign(residualOfY);
          delta_ = deltaY;
        } // Newton Iteration
        jInv_.unbind();
        if( verbose_ )
          std::cerr << "Newton iteration " << iterations_ << ": |residual| = " << delta_ << std::endl;

      }   // operator ()

    };   //class MyNewtonInverseOperator
  }   // namespace Fem

}   // namespace Dune



#endif // DUNE_FEM_MODIFIED_NEWTON_HH
