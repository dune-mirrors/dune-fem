#ifndef DUNE_FEM_SOLVERPARAMETER_HH
#define DUNE_FEM_SOLVERPARAMETER_HH

#include <dune/fem/io/parameter.hh>
#include <dune/common/std/optional.hh>

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
        else if(! verbose_ )
        {
          verbose_ = parameter_.getValue< bool >( keyPrefix_ + "verbose", false );
        }
        return verbose_.value();
      }

      virtual void setVerbose( const bool verb )
      {
        if( other_ )
          other_->setVerbose( verb );
        else
          verbose_ = verb;
      }

      virtual int errorMeasure() const
      {
        if( other_ )
          return other_->errorMeasure();
        else if( !errorMeasure_)
        {
          const std::string errorTypeTable[] =
          { "absolute", "relative", "residualreduction" };
          const int errorType = parameter_.getEnum( keyPrefix_ + "errormeasure", errorTypeTable, 0 );
          errorMeasure_ = errorType ;
        }
        return errorMeasure_.value();
      }

      virtual void setErrorMeasure( const int errorType )
      {
        if( other_ )
          return other_->setErrorMeasure( errorType );
        else
          errorMeasure_ = errorType;
      }

      virtual double linAbsTol ( double eps = 1e-8 )  const
      {
        if( other_ )
          return other_->linAbsTol();
        else if( ! linAbsTol_ )
          linAbsTol_ = parameter_.getValue< double >(keyPrefix_ +  "linabstol", eps );
        return linAbsTol_.value();
      }

      virtual void setLinAbsTol ( const double eps )
      {
        if( other_ )
          return other_->setLinAbsTol( eps );
        else
          linAbsTol_ = eps;
      }

      virtual double linReduction ( double eps = 1e-2 ) const
      {
        if( other_ )
          return other_->linReduction();
        else if ( !linReduction_ )
          linReduction_ = parameter_.getValue< double >( keyPrefix_ + "linreduction", eps );
        return linReduction_.value();
      }

      virtual void setLinReduction ( const double eps )
      {
        if( other_ )
          return other_->setLinReduction( eps );
        else
          linReduction_ = eps;
      }

      virtual int maxLinearIterations () const
      {
        if( other_ )
          return other_->maxLinearIterations();
        else if (! maxLinearIterations_ )
          maxLinearIterations_ = parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
        return maxLinearIterations_.value();
      }

      virtual void  setMaxLinearIterations ( const int maxIter )
      {
        if( other_ )
          return other_->setMaxLinearIterations( maxIter );
        else
          maxLinearIterations_ = maxIter;
      }

      virtual int krylovMethod() const
      {
        if( other_ )
          return other_->krylovMethod();
        else if (! krylovMethod_ )
        {
          const std::string krylovMethodTable[] =
          { "cg", "bicgstab", "gmres", "minres", "gradient", "loop"  };
          int methodType = gmres;
          if( parameter_.exists( keyPrefix_ + "krylovmethod" ) )
            methodType = parameter_.getEnum( keyPrefix_ + "krylovmethod", krylovMethodTable, gmres );
          else
            methodType = parameter_.getEnum( "krylovmethod", krylovMethodTable, gmres );
          krylovMethod_ =  methodType;
        }
        return krylovMethod_.value();
      }

      virtual void setKrylovMethod( const int method )
      {
        if( other_ )
          return other_->setKrylovMethod( method);
        else
          krylovMethod_ = method;
      }

      virtual int gmresRestart() const
      {
        if( other_ )
          return other_->gmresRestart();
        else if (! gmresRestart_ )
        {
          int defaultRestart = 20;
          gmresRestart_ =  parameter_.getValue< int >( keyPrefix_ + "gmres.restart", defaultRestart );
        }
        return gmresRestart_.value();
      }

      virtual void setGmresRestart( const int restart )
      {
        if( other_ )
          return other_->setGmresRestart( restart );
        else
          gmresRestart_ = restart;
      }

      virtual int preconditionMethod () const
      {
        if( other_ )
          return other_ -> preconditionMethod();
        else if(!preconditionMethod_)
        {
          static const std::string preConTable[]
            = { "none", "ssor", "sor", "ilu", "gauss-seidel", "jacobi", "amg-ilu", "amg-jacobi", "ildl" };
          preconditionMethod_ =  parameter_.getEnum(  keyPrefix_ + "preconditioning.method", preConTable, 0 );
        }
        return preconditionMethod_.value();
      }

      virtual void setPreconditionMethod ( const int precMethod )
      {
        if( other_ )
          return other_-> setPreconditionMethod( precMethod );
        else
          preconditionMethod_ = precMethod;
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
        else if (! relaxation_ )
          relaxation_ = parameter_.getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
        return relaxation_.value();
      }

      virtual void setRelaxation ( const double relaxation )
      {
        if( other_ )
          return other_->setRelaxation( relaxation );
        else
          relaxation_ = relaxation;
      }


      virtual int preconditionerIteration () const
      {
        if( other_ )
          return other_->preconditionerIteration();
        else if (! preconditionerIteration_ )
        {
          // TODO: add also check for level
          preconditionerIteration_ =  parameter_.getValue< int >( keyPrefix_ + "preconditioning.iterations", 0 );
        }
        return preconditionerIteration_.value();
      }

      virtual void setPreconditionerIteration ( const int precIter)
      {
        if( other_ )
          return other_->setPreconditionerIteration( precIter );
        else
          preconditionerIteration_ = precIter;
      }


      virtual SolverParameter* clone () const { return new SolverParameter(); }

      //deprecated methods
      [[deprecated]]
      virtual double linAbsTolParameter ()  const
      {
        if( other_ )
          return other_->linAbsTolParameter();
        else
          return parameter_.getValue< double >(keyPrefix_ +  "linabstol", 1e-8 );
      }

      [[deprecated]]
      virtual double linReductionParameter () const
      {
        if( other_ )
          return other_->linReductionParameter();
        else
          return parameter_.getValue< double >( keyPrefix_ + "linreduction", 1e-2 );
      }

      [[deprecated]]
      virtual int maxLinearIterationsParameter () const
      {
        if( other_ )
          return other_->maxLinearIterationsParameter();
        else
          return parameter_.getValue< int >( keyPrefix_ + "maxlineariterations", std::numeric_limits< int >::max() );
      }

    private:
      mutable Std::optional<bool> verbose_;
      mutable Std::optional<int> errorMeasure_;
      mutable Std::optional<double> linAbsTol_;
      mutable Std::optional<double> linReduction_;
      mutable Std::optional<int> maxLinearIterations_;
      mutable Std::optional<int> krylovMethod_;
      mutable Std::optional<int> gmresRestart_;
      mutable Std::optional<int> preconditionMethod_;
      mutable Std::optional<double> relaxation_;
      mutable Std::optional<int> preconditionerIteration_;
    };

  }
}

#endif // #ifndef DUNE_FEM_SOLVERPARAMETER_HH
