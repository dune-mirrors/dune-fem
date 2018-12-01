#ifndef DUNE_FEM_SOLVERPARAMETER_HH
#define DUNE_FEM_SOLVERPARAMETER_HH

#include <dune/fem/io/parameter.hh>

namespace Dune
{

  namespace Fem
  {
    struct SolverParameter
#ifndef DOXYGEN
    : public LocalParameter< SolverParameter, SolverParameter >
#endif
    {
    protected:
      // key prefix, default is fem.solver (can be overloaded by user)
      const std::string keyPrefix_;

      ParameterReader parameter_;

      std::shared_ptr< SolverParameter > other_;

    public:
      // identifier for Fem, ISTL and Petsc solvers
      static const int cg       = 0 ; // CG
      static const int bicgstab = 1 ; // BiCGStab
      static const int gmres    = 2 ; // GMRES
      static const int minres   = 3 ; // MinRes
      static const int gradient = 4 ; // GradientSolver
      static const int loop     = 5 ; // LoopSolver
      static const int superlu  = 6 ; // SuperLUSolver

      static const int none         = 0 ; // no preconditioner
      static const int ssor         = 1 ; // SSOR preconditioner
      static const int sor          = 2 ; // SOR preconditioner
      static const int ilu          = 3 ; // ILU preconditioner
      static const int gauss_seidel = 4 ; // Gauss-Seidel preconditioner
      static const int jacobi       = 5 ; // Jacobi preconditioner
      static const int amg_ilu      = 6 ; // AMG with ILU-0 smoother (deprecated)
      static const int amg_jacobi   = 7 ; // AMG with Jacobi smoother
      static const int ildl         = 8 ; // ILDL from istl

      explicit SolverParameter ( const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( "fem.solver." ), parameter_( parameter ), other_()
      {
      }

      //! constructor passing other implementation as interface which is needed
      //to overload internal default implementations
      SolverParameter ( const SolverParameter* other )
        : keyPrefix_( other->keyPrefix_ ),
          parameter_( other->parameter() ),
          other_()
      {
        // if other is a derived type
        // then store a copy of that object
        if( typeid(SolverParameter) != typeid(*other) )
        {
          //std::cout << typeid(*other).name() << std::endl;
          other_.reset( other->clone() );
        }
      }

      explicit SolverParameter ( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( keyPrefix ), parameter_( parameter ), other_()
      {}

      const ParameterReader& parameter() const { return parameter_; }

      virtual bool verbose() const
      {
        if( other_ )
          return other_->verbose();
        else
          return parameter_.getValue< bool >( keyPrefix_ + "verbose", false );
      }

      virtual int errorMeasure() const
      {
        if( other_ )
          return other_->errorMeasure();
        else
        {
          const std::string errorTypeTable[] =
            { "absolute", "relative", "residualreduction" };
          const int errorType = parameter_.getEnum( keyPrefix_ + "errormeasure", errorTypeTable, 0 );
          return errorType ;
        }
      }

      virtual double linAbsTolParameter ()  const
      {
        if( other_ )
          return other_->linAbsTolParameter();
        else
          return parameter_.getValue< double >(keyPrefix_ +  "linabstol", 1e-8 );
      }

      virtual double linReductionParameter () const
      {
        if( other_ )
          return other_->linReductionParameter();
        else
          return parameter_.getValue< double >( keyPrefix_ + "linreduction", 1e-2 );
      }

      virtual int maxLinearIterationsParameter () const
      {
        if( other_ )
          return other_->maxLinearIterationsParameter();
        else
          return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
      }

      virtual int krylovMethod() const
      {
        if( other_ )
          return other_->krylovMethod();
        else
        {
          const std::string krylovMethodTable[] =
            { "cg", "bicgstab", "gmres", "minres", "gradient", "loop"  };
          int methodType = gmres;
          if( parameter_.exists( keyPrefix_ + "krylovmethod" ) )
            methodType = parameter_.getEnum( keyPrefix_ + "krylovmethod", krylovMethodTable, gmres );
          else
            methodType = parameter_.getEnum( "krylovmethod", krylovMethodTable, gmres );
          return methodType;
        }
      }

      virtual int gmresRestart() const
      {
        if( other_ )
          return other_->gmresRestart();
        else
        {
          int defaultRestart = 20;
          return parameter_.getValue< int >( keyPrefix_ + "gmres.restart", defaultRestart );
        }
      }

      virtual int preconditionMethod () const
      {
        static const std::string preConTable[]
          = { "none", "ssor", "sor", "ilu", "gauss-seidel", "jacobi", "amg-ilu", "amg-jacobi", "ildl" };
        return parameter_.getEnum(  keyPrefix_ + "preconditioning.method", preConTable, 0 );
      }

      virtual std::string preconditionName() const
      {
        static const std::string preConTable[]
          = { "None", "SSOR", "SOR", "ILU", "Gauss-Seidel", "Jacobi", "AMG-ILU", "AMG-Jacobi", "ILDL" };
        const int precond = preconditionMethod();
        std::stringstream tmp;
        tmp << preConTable[precond];

        if( precond != 3 )
          tmp << " n=" << preconditionerIteration();
        tmp << " relax=" << relaxation();
        return tmp.str();
      }

      virtual double relaxation () const
      {
        if( other_ )
          return other_->relaxation();
        else
          return parameter_.getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
      }

      virtual int preconditionerIteration () const
      {
        if( other_ )
          return other_->preconditionerIteration();
        else
        {
          // TODO: add also check for level
          return parameter_.getValue< int >( keyPrefix_ + "preconditioning.iterations", 0 );
        }
      }

      virtual SolverParameter* clone () const { return new SolverParameter(); }
    };

  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
