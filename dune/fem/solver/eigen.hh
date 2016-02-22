#ifndef DUNE_FEM_SOLVER_EIGEN_HH
#define DUNE_FEM_SOLVER_EIGEN_HH

#ifdef HAVE_EIGEN

#include <Eigen/IterativeLinearSolvers>
#include <dune/fem/operator/linear/eigenoperator.hh>
#include <dune/fem/operator/common/operator.hh>

namespace Dune
{

  namespace Fem
  {
    template< class DiscreteFunction >
    class EigenCGInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;

    private:
      typedef typename OperatorType::MatrixType::MatrixStorageType Matrix;

    public:
      EigenCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : op_(op)
        , solver_()
      {
      }

      EigenCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : op_(op)
        , solver_()
      {
        solver_.setTolerance(absLimit);
        solver_.analyzePattern( op_.matrix().data() );
        solver_.factorize( op_.matrix().data() );
      }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        std::cout << "starting to solve using eigen" << std::endl;
        for (int i=0;i<100;++i)
        {
          std::cout << "starting to solve" << std::endl;
          w.dofVector().array().coefficients()
             = solver_.solve( u.dofVector().array().coefficients() );
        }
      }

      unsigned int iterations () const { return solver_.iterations(); }

    protected:
      const OperatorType &op_;
      Eigen::ConjugateGradient<Matrix> solver_;
    };
  } // namespace Fem

} // namespace Dune

#endif
#endif // #ifndef DUNE_FEM_SOLVER_VIENNACL_HH
