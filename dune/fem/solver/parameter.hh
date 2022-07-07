#ifndef DUNE_FEM_SOLVERPARAMETER_HH
#define DUNE_FEM_SOLVERPARAMETER_HH

#include <dune/fem/io/parameter.hh>

namespace Dune
{

  namespace Fem
  {
    namespace LinearSolver
    {
      struct ToleranceCriteria {
        static const int absolute = 0;
        static const int relative = 1;
        static const int residualReduction = 2;
      };
    }

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
      // identifier for Fem, ISTL and Petsc solvers
      static const int cg       = 1 ; // CG
      static const int bicgstab = 2 ; // BiCGStab
      static const int gmres    = 3 ; // GMRES
      static const int minres   = 4 ; // MinRes
      static const int gradient = 5 ; // GradientSolver
      static const int loop     = 6 ; // LoopSolver
      static const int superlu  = 7 ; // SuperLUSolver
      static const int bicg     = 8 ; // BiCG
      static const int preonly  = 9 ; // only preconder
      static const std::string solverMethodTable(int i)
      {
        std::string methods[] =
           { "cg", "bicgstab", "gmres", "minres", "gradient", "loop", "superlu", "bicg", "preonly"  };
        return methods[i-1]; // starting to count from 1,...
      }

      static const int none         = 1 ; // no preconditioner
      static const int ssor         = 2 ; // SSOR preconditioner
      static const int sor          = 3 ; // SOR preconditioner
      static const int ilu          = 4 ; // ILU preconditioner
      static const int gauss_seidel = 5 ; // Gauss-Seidel preconditioner
      static const int jacobi       = 6 ; // Jacobi preconditioner
      static const int amg_ilu      = 7 ; // AMG with ILU-0 smoother (deprecated)
      static const int amg_jacobi   = 8 ; // AMG with Jacobi smoother
      static const int ildl         = 9 ; // ILDL from istl
      static const int oas          = 10; // Overlapping Additive Schwarz
      static const int icc          = 11; // Incomplete Cholesky factorization
      static const std::string preconditionMethodTable(int i)
      {
        static const std::string methods[]
          = { "none", "ssor", "sor", "ilu", "gauss-seidel", "jacobi", "amg-ilu", "amg-jacobi", "ildl", "oas", "icc" };
        return methods[i-1]; // starting to count from 1,...
      }

      SolverParameter ( const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( "fem.solver." ), parameter_( parameter )
      {}

      explicit SolverParameter ( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( keyPrefix ), parameter_( parameter )
      {}

      const std::string& keyPrefix() const { return keyPrefix_; }

      const ParameterReader& parameter() const { return parameter_; }

      virtual void reset()
      {
        verbose_       = -1;
        absoluteTol_   = -1;
        reductionTol_  = -1;
        maxIterations_ = -1;
      }

      virtual bool verbose() const
      {
        if( verbose_ < 0 )
        {
          verbose_ = parameter_.getValue< bool >( keyPrefix_ + "verbose", false ) ? 1 : 0;
        }
        return bool(verbose_);
      }

      virtual void setVerbose( const bool verb )
      {
        verbose_ = verb ? 1 : 0;
      }

      virtual int errorMeasure() const
      {
        const std::string errorTypeTable[] =
        { "absolute", "relative", "residualreduction" };
        const int errorType = parameter_.getEnum( keyPrefix_ + "errormeasure", errorTypeTable, 0 );
        return errorType ;
      }

      virtual double tolerance (  ) const
      {
        if( tolerance_ < 0 )
        {
          double defaultTol = 1e-8;
          if(parameter_.exists(keyPrefix_ + "tolerance"))
            tolerance_ = parameter_.getValue< double >( keyPrefix_ + "tolerance", defaultTol );
          else
          {
            if( parameter_.exists(keyPrefix_ + "absolutetol") ||
                parameter_.exists(keyPrefix_ + "reductiontol") )
              // new parameter not used but old parameters exist
            {
              int measure = errorMeasure();
              if (measure == 0)
                tolerance_ = absoluteTol__();
              else
                tolerance_ = reductionTol__();
            }
            else tolerance_ = defaultTol; // no parameter set
          }
        }
        return tolerance_;
      }

      virtual void setTolerance ( const double eps )
      {
        assert( eps >= 0.0 );
        tolerance_ = eps;
      }

      virtual int maxIterations () const
      {
        if( maxIterations_ < 0 )
        {
          if(parameter_.exists(keyPrefix_ + "maxlineariterations"))
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ +"maxlineariterations is deprecated. Please use " + keyPrefix_ + "maxiterations instead." << std::endl;
            maxIterations_ =  parameter_.getValue< double >(keyPrefix_ + "maxlineariterations");
          }
          else
            maxIterations_ =  parameter_.getValue< int >( keyPrefix_ + "maxiterations", std::numeric_limits< int >::max() );
        }
        return maxIterations_;
      }

      virtual void  setMaxIterations ( const int maxIter )
      {
        assert( maxIter >= 0 );
        maxIterations_ = maxIter;
      }

      virtual int solverMethod(
            const std::vector<int> standardMethods,
            const std::vector<std::string> &additionalMethods = {},
            int defaultMethod = 0   // this is the first method passed in
          ) const
      {
        std::vector<std::string> methodTable(standardMethods.size()+additionalMethods.size());
        for (std::size_t i=0;i<standardMethods.size();++i)
          methodTable[i] = solverMethodTable(standardMethods[i]);
        for (std::size_t i=0;i<additionalMethods.size();++i)
          methodTable[standardMethods.size()+i] = additionalMethods[i];
        std::size_t method;
        if( parameter_.exists( keyPrefix_ + "method" ) ||
           !parameter_.exists( "method" ) )
          method = parameter_.getEnum( keyPrefix_ + "method", methodTable, defaultMethod );
        else
        {
          method = parameter_.getEnum( "krylovmethod", methodTable, defaultMethod );
          std::cout << "WARNING: using old parameter name 'krylovmethod' "
                    << "please switch to '" << keyPrefix_ << "method'\n";
        }
        if (method < standardMethods.size())
          return standardMethods[method];
        else
          return -(method-standardMethods.size()); // return in [ 0,-additionalMethods.size() )
      }

      virtual int gmresRestart() const
      {
        int defaultRestart = 20;
        return parameter_.getValue< int >( keyPrefix_ + "gmres.restart", defaultRestart );
      }

      virtual int preconditionMethod(
            const std::vector<int> standardMethods,
            const std::vector<std::string> &additionalMethods = {},
            int defaultMethod = 0   // this is the first method passed in
          ) const
      {
        std::vector<std::string> methodTable(standardMethods.size()+additionalMethods.size());
        for (std::size_t i=0;i<standardMethods.size();++i)
          methodTable[i] = preconditionMethodTable(standardMethods[i]);
        for (std::size_t i=0;i<additionalMethods.size();++i)
          methodTable[standardMethods.size()+i] = additionalMethods[i];
        std::size_t method = parameter_.getEnum( keyPrefix_ + "preconditioning.method", methodTable, defaultMethod );
        if (method < standardMethods.size())
          return standardMethods[method];
        else
          return -(method-standardMethods.size()); // return in [ 0,-additionalMethods.size() )
      }

      virtual double relaxation () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
      }

      virtual int preconditionerIteration () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "preconditioning.iterations", 1 );
      }

      virtual int preconditionerLevel () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "preconditioning.level", 0 );
      }

      virtual bool threading () const
      {
        return parameter_.getValue< bool >( keyPrefix_ + "threading", true );
      }

     private:
      virtual double absoluteTol__ ( )  const
      {
        if( absoluteTol_ < 0 )
        {
          if(parameter_.exists(keyPrefix_ + "linabstol"))
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ + "linabstol is deprecated. Please use " + keyPrefix_ + "absolutetol instead." << std::endl;
            absoluteTol_ =  parameter_.getValue< double >(keyPrefix_ + "linabstol");
          }
          else
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ + "absolutetol is deprecated. Please use " + keyPrefix_ + "tolerance instead." << std::endl;
            absoluteTol_ =  parameter_.getValue< double >(keyPrefix_ +  "absolutetol", 1e-8 );
          }
        }
        return absoluteTol_;
      }

      virtual double reductionTol__ (  ) const
      {
        if( reductionTol_ < 0 )
        {
          if(parameter_.exists(keyPrefix_ + "linreduction"))
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ +"linreduction is deprecated. Please use " + keyPrefix_ + "reductiontol instead." << std::endl;
            reductionTol_ =  parameter_.getValue< double >(keyPrefix_ + "linreduction");
          }
          else
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ +"reductiontol is deprecated. Please use " + keyPrefix_ + "tolerance instead." << std::endl;
            reductionTol_ = parameter_.getValue< double >( keyPrefix_ + "reductiontol", 1e-8 );
          }
        }
        return reductionTol_;
      }

      mutable int    verbose_       = -1;
      mutable int    maxIterations_ = -1;
      mutable double absoluteTol_   = -1.;
      mutable double reductionTol_  = -1.;
      mutable double tolerance_     = -1.;
    };
  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
