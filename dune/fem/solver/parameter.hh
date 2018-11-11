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

    public:
      // identifier for Fem and ISTL solvers
      static const int cg       = 0 ; // CG
      static const int bicgstab = 1 ; // BiCGStab
      static const int gmres    = 2 ; // GMRES
      static const int minres   = 3 ; // MinRes
      static const int gradient = 4 ; // GradientSolver
      static const int loop     = 5 ; // LoopSolver
      static const int superlu  = 6 ; // SuperLUSolver

      explicit SolverParameter ( const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( "fem.solver." ), parameter_( parameter )
      {
      }

      explicit SolverParameter ( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( keyPrefix ), parameter_( parameter )
      {}

      const ParameterReader& parameter() const { return parameter_; }

      virtual bool verbose() const
      {
        return parameter_.getValue< bool >( keyPrefix_ + "verbose", false );
      }

      virtual int errorMeasure() const
      {
        const std::string errorTypeTable[] =
          { "absolute", "relative", "residualreduction" };
        const int errorType = parameter_.getEnum( keyPrefix_ + "errormeasure", errorTypeTable, 0 );
        return errorType ;
      }

      virtual double linAbsTolParameter ()  const
      {
        return parameter_.getValue< double >(keyPrefix_ +  "linabstol", 1e-8 );
      }

      virtual double linReductionParameter () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "linreduction", 1e-2 );
      }

      virtual int maxLinearIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
      }

      virtual int krylovMethod() const
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

      virtual int gmresRestart() const
      {
        int defaultRestart = 20;
        return parameter_.getValue< int >( keyPrefix_ + "gmres.restart", defaultRestart );
      }

    };

  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
