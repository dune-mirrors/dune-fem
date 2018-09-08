#ifndef DUNE_FEM_AMGXSOLVER_HH
#define DUNE_FEM_AMGXSOLVER_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#if HAVE_PETSC_AMGX
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/solver/parameter.hh>

#include <AmgXSolver.hpp>

namespace Dune
{

  namespace Fem
  {

    //=====================================================================
    // Implementation of AMGX solver wrapper using PETSc matrix
    //=====================================================================

    /** @ingroup OEMSolver
        @{
    **/

    // AMGXSolver
    // --------------

    /** \brief AMGX solver context for PETSc Mat and PETSc Vec */
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    class AMGXInverseOperator
      : public Operator< DF, DF >
    {
    public:
      typedef DF DiscreteFunctionType;
      typedef DiscreteFunctionType  DestinationType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef PetscDiscreteFunction< DiscreteFunctionSpaceType > PetscDiscreteFunctionType;

      typedef PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType >  AssembledOperatorType;

      typedef Op OperatorType;

      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolut limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note AMGX solvers uses the relative reduction.
       */

      AMGXInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             int maxIter,
                             bool verbose,
                             const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        parameter_( parameter )
      {
        bind( op );
      }

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      AMGXInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             int maxIter,
                             const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {
        bind( op );
      }

      AMGXInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( std::numeric_limits< int >::max()),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {
        bind( op );
      }

      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolut limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note PETSc Krylov solvers uses the relative reduction.
       */

      AMGXInverseOperator ( double reduction, double absLimit, int maxIter, bool verbose, const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        parameter_( parameter )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      AMGXInverseOperator ( double reduction, double absLimit, int maxIter, const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {}

      AMGXInverseOperator ( double reduction, double absLimit, const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( std::numeric_limits< int >::max()),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {}

      AMGXInverseOperator ( const OperatorType &op,
          double reduction, double absLimit,
          const SolverParameter &parameter )
      : AMGXInverseOperator( reduction, absLimit, parameter.parameter() )
      {}

      AMGXInverseOperator ( double reduction, double absLimit,
          const SolverParameter &parameter )
      : AMGXInverseOperator( reduction, absLimit, parameter.parameter() )
      {}

      void bind ( const OperatorType &op )
      {
        op_ = &op;
        matrixOp_ = dynamic_cast<const AssembledOperatorType*> ( &op );
        initialize( parameter_ );
      }

      void unbind ()
      {
        op_ = nullptr; matrixOp_ = nullptr;
        amgXSolver_.finalize();
      }

      void initialize ( const ParameterReader &parameter )
      {
        if( !matrixOp_ )
          DUNE_THROW(NotImplemented,"AMGX solver with matrix free implementations is not supported!");

        std::string mode = "dDDI";
        std::string solverconfig = "./amgxconfig";

        // matrixOp_->space().gridPart.comm().communicator()
        amgXSolver_.initialize(PETSC_COMM_WORLD, mode, solverconfig);

        // attach Matrix to linear solver context
        Mat& A = const_cast<Mat &> (matrixOp_->petscMatrix());

        // set matrix
        amgXSolver_.setA( A );
      }

      void setMaxIterations ( std::size_t maxIter )
      {
        maxIter_ = maxIter;
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {}

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: " << solverName_ << " eps = " << reduction_ ;
        out  << "\\\\ \n";
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      template <class DiscreteFunction>
      void apply( const DiscreteFunction& arg, DiscreteFunction& dest ) const
      {
        // copy discrete functions
        PetscDiscreteFunctionType Arg("PetscSolver::arg", arg.space() );
        Arg.assign( arg );

        // also copy initial destination in case this is used a solver init value
        PetscDiscreteFunctionType Dest("PetscSolver::dest", dest.space() );
        Dest.assign( dest );

        apply( Arg, Dest );
        // copy destination back
        dest.assign( Dest );
      }

      void apply( const PetscDiscreteFunctionType& arg, PetscDiscreteFunctionType& dest ) const
      {
        // need to have a 'distributed' destination vector for continuous spaces
        if( dest.space().continuous() )
          dest.dofVector().clearGhost();

        // call PETSc solvers, dest = x, arg = rhs
        Vec& x = *dest.petscVec();
        Vec& rhs = *(const_cast< PetscDiscreteFunctionType& > (arg).petscVec());
        amgXSolver_.solve( x, rhs );

        // a continuous solution is 'distributed' so need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

        // get number of iterations
        amgXSolver_.getIters( iterations_ );
      }

      // return number of iterations
      int iterations() const
      {
        return iterations_;
      }

      //! return accumulated communication time
      double averageCommTime() const
      {
        return -1;
      }

      /** \brief solve the system
          \param[in] arg right hand side
          \param[out] dest solution
      */
      void operator() ( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

      template <class DiscreteFunction>
      void operator() ( const DiscreteFunction& arg, DiscreteFunction& dest ) const
      {
        apply(arg,dest);
      }

    protected:
      const OperatorType  *op_ = nullptr; // linear operator
      const AssembledOperatorType* matrixOp_ = nullptr; // assembled operator

      mutable AmgXSolver amgXSolver_;

      double reduction_;
      double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_ = 0;
      ParameterReader parameter_;
      std::string solverName_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCSOLVER_HH
