#ifndef DUNE_FEM_SOLVER_EIGEN_HH
#define DUNE_FEM_SOLVER_EIGEN_HH

#ifdef HAVE_EIGEN

#include <Eigen/IterativeLinearSolvers>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>
#include <dune/fem/solver/parameter.hh>

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
      EigenInverseOperator ( double redEps, double absLimit, unsigned int maxIter, bool verbose,
                             const ParameterReader &parameter = Parameter::container() )
        : solver_(std::make_unique<EigenOp>()), absLimit_( absLimit )
      {}

      EigenInverseOperator ( double redEps, double absLimit, unsigned int maxIter = std::numeric_limits<int>::max(),
                             const ParameterReader &parameter = Parameter::container() )
        : solver_(std::make_unique<EigenOp>()), absLimit_( absLimit )
      {}

      EigenInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose,
                             const ParameterReader &parameter = Parameter::container() )
        : solver_(std::make_unique<EigenOp>()), absLimit_( absLimit )
      {
        bind( op );
      }

      EigenInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter,
                             const ParameterReader &parameter = Parameter::container() )
        : solver_(std::make_unique<EigenOp>()), absLimit_( absLimit )
      {
        bind( op );
      }

      EigenInverseOperator ( const OperatorType &op,
          double reduction, double absLimit,
          const SolverParameter &parameter )
      : EigenInverseOperator( reduction, absLimit, std::numeric_limits<int>::max(), parameter.parameter() )
      {}
      EigenInverseOperator ( double reduction, double absLimit,
          const SolverParameter &parameter )
      : EigenInverseOperator( reduction, absLimit, std::numeric_limits<int>::max(), parameter.parameter() )
      {}

      void bind ( const OperatorType &op )
      {
        op_ = &op;
        setup();
      }
      void unbind() { op_ = nullptr; }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        if ( op_ )
          w.dofVector().array().coefficients() = solver_->solve( u.dofVector().array().coefficients() );
      }

      unsigned int iterations () const { return solver_->iterations(); }
      void setMaxIterations ( unsigned int ) {}

    protected:
      void setup ()
      {
        assert( op_ );
        solver_->setTolerance( absLimit_ );
        solver_->analyzePattern( op_->matrix().data() );
        solver_->factorize( op_->matrix().data() );
      }

      const OperatorType *op_ = nullptr;
      std::unique_ptr<EigenOp> solver_;
      double absLimit_;
    };

    template< class DiscreteFunction>
    using EigenCGInverseOperator = EigenInverseOperator<DiscreteFunction,
            Eigen::ConjugateGradient<typename EigenInverseOperatorTraits<DiscreteFunction>::Matrix> >;

    template< class DiscreteFunction>
    using EigenBiCGStabInverseOperator = EigenInverseOperator<DiscreteFunction,
            Eigen::BiCGSTAB<typename EigenInverseOperatorTraits<DiscreteFunction>::Matrix> >;

  } // namespace Fem
} // namespace Dune

#endif
#endif // #ifndef DUNE_FEM_SOLVER_VIENNACL_HH
