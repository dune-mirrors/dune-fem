#ifndef DUNE_FEM_SOLVER_EIGEN_HH
#define DUNE_FEM_SOLVER_EIGEN_HH

#ifdef HAVE_EIGEN

#include <Eigen/IterativeLinearSolvers>

#include <dune/fem/operator/linear/eigenoperator.hh>

namespace Dune
{

  namespace Fem
  {
    template <class DiscreteFunction>
    struct EigenInverseOperatorTraits
    {
      typedef Fem::EigenLinearOperator< DiscreteFunction, DiscreteFunction> OperatorType;
      typedef typename OperatorType::MatrixType::MatrixStorageType Matrix;
    };

    template< class DiscreteFunction, class EigenOp >
    class EigenInverseOperator
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
      template <class P>
      EigenInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose,
          const P &p=P())
        : op_(op)
        , solver_()
      {
        solver_.setTolerance(absLimit);
        solver_.analyzePattern( op_.matrix().data() );
        solver_.factorize( op_.matrix().data() );
      }

      template <class P>
      EigenInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter,
          const P &p=P())
        : op_(op)
        , solver_()
      {
        solver_.setTolerance(absLimit);
        solver_.analyzePattern( op_.matrix().data() );
        solver_.factorize( op_.matrix().data() );
      }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        w.dofVector().array().coefficients()
           = solver_.solve( u.dofVector().array().coefficients() );
      }

      unsigned int iterations () const { return solver_.iterations(); }

    protected:
      const OperatorType &op_;
      EigenOp solver_;
    };

    template< class DiscreteFunction>
    class EigenCGInverseOperator
      : public EigenInverseOperator<DiscreteFunction,
            Eigen::ConjugateGradient<typename EigenInverseOperatorTraits<DiscreteFunction>::Matrix>
          >
    {
      typedef EigenInverseOperatorTraits<DiscreteFunction> Traits;
      typedef EigenInverseOperator<DiscreteFunction,
            Eigen::ConjugateGradient<typename Traits::Matrix> > BaseType;
      public:
      typedef typename Traits::OperatorType OperatorType;

      template <class P>
      EigenCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose,
          const P& p=P())
        : BaseType(op,redEps,absLimit,maxIter,verbose,p) {}
      template <class P>
      EigenCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter,
          const P& p=P())
        : BaseType(op,redEps,absLimit,maxIter) {}
      EigenCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : BaseType(op,redEps,absLimit,maxIter) {}
    };
    template< class DiscreteFunction>
    class EigenBiCGStabInverseOperator
      : public EigenInverseOperator<DiscreteFunction,
            Eigen::BiCGSTAB<typename EigenInverseOperatorTraits<DiscreteFunction>::Matrix>
          >
    {
      typedef EigenInverseOperatorTraits<DiscreteFunction> Traits;
      typedef EigenInverseOperator<DiscreteFunction,
            Eigen::BiCGSTAB<typename Traits::Matrix> > BaseType;
      public:
      typedef typename Traits::OperatorType OperatorType;

      template <class P>
      EigenBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, int maxIter, bool verbose,
          const P& p=P())
        : BaseType(op,redEps,absLimit,maxIter,verbose,p) {}
      template <class P>
      EigenBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter,
          const P& p=P())
        : BaseType(op,redEps,absLimit,maxIter) {}
      EigenBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : BaseType(op,redEps,absLimit,maxIter) {}
    };

  } // namespace Fem
} // namespace Dune

#endif
#endif // #ifndef DUNE_FEM_SOLVER_VIENNACL_HH
