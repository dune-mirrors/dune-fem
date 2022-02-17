#ifndef DUNE_FEM_AMGXSOLVER_HH
#define DUNE_FEM_AMGXSOLVER_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/solver/inverseoperatorinterface.hh>
#include <dune/fem/solver/parameter.hh>
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/function/petscdiscretefunction.hh>

#if HAVE_AMGXSOLVER
// AMGX solver wrapper based on Petsc data structures
#include <AmgXSolver.hpp>
#endif

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

    struct AMGXSolverParameter : public LocalParameter< SolverParameter, AMGXSolverParameter >
    {
      typedef LocalParameter< SolverParameter, AMGXSolverParameter >  BaseType;

    public:
      using BaseType :: parameter;
      using BaseType :: keyPrefix;

      AMGXSolverParameter( const ParameterReader& parameter = Parameter::container() )
        : BaseType( parameter )
      {}


      AMGXSolverParameter( const std::string &keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : BaseType( keyPrefix, parameter )
      {}

      AMGXSolverParameter( const SolverParameter& sp )
        : AMGXSolverParameter( sp.keyPrefix(), sp.parameter() )
      {}

      virtual std::string solvermode () const
      {
        const std::string modes [] = { "dDDI" , "dDFI", "dFFI", "hDDI", "hDFI", "hFFI" };
        int mode = parameter().getEnum(keyPrefix() + "amgx.mode", modes, 0 );
        return modes[ mode ];
      }

      virtual std::string solverconfig () const
      {
        return parameter().template getValue< std::string >( keyPrefix() + "amgx.config", "amgxconfig.json");
      }
    };



    // AMGXSolver
    // --------------

    template< class DiscreteFunction >
    class AMGXInverseOperator;

    template< class DiscreteFunction >
    struct AMGXInverseOperatorTraits
    {
      typedef DiscreteFunction    DiscreteFunctionType;
      typedef PetscDiscreteFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > SolverDiscreteFunctionType;

      typedef Dune::Fem::Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef OperatorType  PreconditionerType;

#if HAVE_AMGXSOLVER
      typedef Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction > AssembledOperatorType;
#else
      typedef OperatorType  AssembledOperatorType;
#endif

      typedef AMGXInverseOperator< DiscreteFunction >  InverseOperatorType;

      typedef AMGXSolverParameter SolverParameterType;
    };




    /** \brief AMGX solver context for PETSc Mat and PETSc Vec */
    template< class DF >
    class AMGXInverseOperator
      : public InverseOperatorInterface< AMGXInverseOperatorTraits< DF > >
    {
      typedef AMGXInverseOperatorTraits< DF >  Traits;
      typedef InverseOperatorInterface< Traits >  BaseType;
      friend class InverseOperatorInterface< Traits >;
    public:
      using BaseType :: parameter;

      /** \brief this solver does not offer to set preconditioning option */
      static const bool preconditioningAvailable = false;

      typedef typename BaseType::SolverDiscreteFunctionType  SolverDiscreteFunctionType;
      typedef typename BaseType::OperatorType                OperatorType;
      typedef typename BaseType::PreconditionerType          PreconditionerType;
      typedef typename BaseType::AssembledOperatorType       AssembledOperatorType;

      /** \brief constructor
       *
       *  \param[in] parameter parameter for the solver
       */
      AMGXInverseOperator ( const AMGXSolverParameter& parameter = AMGXSolverParameter() )
        : BaseType( parameter )
      {
      }

      AMGXInverseOperator ( const OperatorType &op,
                            const AMGXSolverParameter & parameter = AMGXSolverParameter() )
        : AMGXInverseOperator ( parameter )
      {
        bind( op );
      }

      AMGXInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                            const AMGXSolverParameter & parameter = AMGXSolverParameter() )
        : AMGXInverseOperator( parameter )
      {
        bind( op, preconditioner );
      }

      AMGXInverseOperator ( const AMGXInverseOperator& other )
        : AMGXInverseOperator( other.parameter() )
      {
        if( other.operator_ )
          bind( *(other.operator_) );
      }

      void bind( const OperatorType& op )
      {
        BaseType::bind( op );
        init( parameter() );
      }

      void unbind()
      {
#if HAVE_AMGXSOLVER
        amgXSolver_->finalize();
        amgXSolver_.reset();
#endif
        BaseType :: unbind();
      }

    protected:
      void init( const AMGXSolverParameter& parameter )
      {
        if( assembledOperator_ )
        {
          std::string mode   = parameter.solvermode();
          std::string config = parameter.solverconfig();
#if HAVE_AMGXSOLVER
          amgXSolver_.reset( new AmgXSolver() );
          amgXSolver_->initialize(PETSC_COMM_WORLD, mode, config );

          // check that PetscMat was assembled not in block mode
          if( assembledOperator_->blockedMode() )
            DUNE_THROW(InvalidStateException, "AMGXInverseOperator only works with PetscLinearOperator in non-blocked mode!");

          // attach Matrix to linear solver context
          Mat& A = const_cast<Mat &> (assembledOperator_->exportMatrix());

          // set matrix
          amgXSolver_->setA( A );
#else
          DUNE_THROW(InvalidStateException,"AMGX solver or PETSc not found during cmake config. Please reconfigure!");
#endif
        }
      }

      int apply( const SolverDiscreteFunctionType& arg, SolverDiscreteFunctionType& dest ) const
      {
        if( !assembledOperator_ )
          DUNE_THROW(NotImplemented,"AMGX solver with matrix free implementations is not supported!");


        int iterations = -1;
#if HAVE_AMGXSOLVER
        assert( amgXSolver_ );

        // need to have a 'distributed' destination vector for continuous spaces
        if( dest.space().continuous() )
          dest.dofVector().clearGhost();

        // call PETSc solvers, dest = x, arg = rhs
        Vec& x = *dest.petscVec();
        Vec& rhs = *(const_cast< SolverDiscreteFunctionType& > (arg).petscVec());
        amgXSolver_->solve( x, rhs );

        // a continuous solution is 'distributed' so need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

        // get number of iterations
        amgXSolver_->getIters( iterations );
#else
        DUNE_THROW(InvalidStateException,"AMGX solver or PETSc not found during cmake config. Please reconfigure!");
#endif
        return iterations;
      }

    protected:
#if HAVE_AMGXSOLVER
      mutable std::unique_ptr< AmgXSolver > amgXSolver_;
#endif
      using BaseType :: assembledOperator_;
      using BaseType :: parameter_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PETSCSOLVER_HH
