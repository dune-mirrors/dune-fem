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
        : keyPrefix_( "fem.solver." ), parameter_( parameter )
      {
      }

      explicit SolverParameter ( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : keyPrefix_( keyPrefix ), parameter_( parameter )
      {}

      const std::string keyPrefix() const { return keyPrefix_; }

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

      virtual void setErrorMeasure( const int errorType )
      {
        const std::string errorTypeTable[] =
        { "absolute", "relative", "residualreduction" };
        Parameter::append( keyPrefix_ + "errormeasure", errorTypeTable[errorType], true );
      }

      virtual double absoluteTol ( )  const
      {
        if( absoluteTol_ < 0 )
        {
          if(parameter_.exists(keyPrefix_ + "linabstol"))
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ + "linabstol is deprecated. Please use " + keyPrefix_ + "absolutetol instead." << std::endl;
            absoluteTol_ =  parameter_.getValue< double >(keyPrefix_ + "linabstol");
          }
          else
            absoluteTol_ =  parameter_.getValue< double >(keyPrefix_ +  "absolutetol", 1e-8 );
        }
        return absoluteTol_;
      }

      virtual void setAbsoluteTol ( const double eps )
      {
        assert( eps >= 0.0 );
        absoluteTol_ = eps;
      }

      virtual double reductionTol (  ) const
      {
        if( reductionTol_ < 0 )
        {
          if(parameter_.exists(keyPrefix_ + "linreduction"))
          {
            std::cout << "WARNING: Parameter " + keyPrefix_ +"linreduction is deprecated. Please use " + keyPrefix_ + "reductiontol instead." << std::endl;
            reductionTol_ =  parameter_.getValue< double >(keyPrefix_ + "linreduction");
          }
          else
            reductionTol_ = parameter_.getValue< double >( keyPrefix_ + "reductiontol", 1e-2 );
        }
        return reductionTol_;
      }

      virtual void setReductionTol ( const double eps )
      {
        assert( eps >= 0.0 );
        reductionTol_ = eps;
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

      virtual void setKrylovMethod( const int method )
      {
        const std::string krylovMethodTable[] =
        { "cg", "bicgstab", "gmres", "minres", "gradient", "loop"  };
        Parameter::append( keyPrefix_ + "krylovmethod", krylovMethodTable[method], true );
      }

      virtual int gmresRestart() const
      {
        int defaultRestart = 20;
        return parameter_.getValue< int >( keyPrefix_ + "gmres.restart", defaultRestart );
      }

      virtual void setGmresRestart( const int restart )
      {
        Parameter::append( keyPrefix_ + "gmres.restart", std::to_string(restart), true );
      }

      virtual int preconditionMethod () const
      {
        static const std::string preConTable[]
          = { "none", "ssor", "sor", "ilu", "gauss-seidel", "jacobi", "amg-ilu", "amg-jacobi", "ildl" };
        return parameter_.getEnum(  keyPrefix_ + "preconditioning.method", preConTable, 0 );
      }

      virtual void setPreconditionMethod ( const int precMethod )
      {
        static const std::string preConTable[]
          = { "none", "ssor", "sor", "ilu", "gauss-seidel", "jacobi", "amg-ilu", "amg-jacobi", "ildl" };
        Parameter::append(  keyPrefix_ + "preconditioning.method", preConTable[precMethod], true );
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
        return parameter_.getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
      }

      virtual void setRelaxation ( const double relaxation )
      {
        Parameter::append( keyPrefix_ + "preconditioning.relaxation", std::to_string(relaxation), true );
      }


      virtual int preconditionerIteration () const
      {
        // TODO: add also check for level
        return parameter_.getValue< int >( keyPrefix_ + "preconditioning.iterations", 0 );
      }

      virtual void setPreconditionerIteration ( const int precIter)
      {
        // TODO: add also check for level
        Parameter::append( keyPrefix_ + "preconditioning.iterations", std::to_string(precIter), true );
      }

      //deprecated methods
      [[deprecated]]
      virtual double linAbsTolParameter ()  const
      {
        return parameter_.getValue< double >(keyPrefix_ +  "linabstol", 1e-8 );
      }

      [[deprecated]]
      virtual double linReductionParameter () const
      {
        return parameter_.getValue< double >( keyPrefix_ + "linreduction", 1e-2 );
      }

      [[deprecated]]
      virtual int maxLinearIterationsParameter () const
      {
        return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
      }

     private:
      mutable int    verbose_       = -1;
      mutable int    maxIterations_ = -1;
      mutable double absoluteTol_   = -1.;
      mutable double reductionTol_  = -1.;
    };

  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
