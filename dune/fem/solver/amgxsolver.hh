#ifndef DUNE_FEM_AMGXSOLVER_HH
#define DUNE_FEM_AMGXSOLVER_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#if HAVE_AMGXSOLVER
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

    struct AMGXSolverParameter : public LocalParameter< SolverParameter, AMGXSolverParameter >
    {
      typedef LocalParameter< SolverParameter, AMGXSolverParameter > BaseType;

      AMGXSolverParameter( const ParameterReader &parameter = Parameter::container() )
        : BaseType( parameter )
      {}

      AMGXSolverParameter( const std::string keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : BaseType( keyPrefix, parameter )
      {}

      AMGXSolverParameter( const SolverParameter& other )
        : BaseType( other.keyPrefix(), other.parameter() )
      {}


      AMGXSolverParameter( const std::string keyPrefix = "fem.solver." )
        : BaseType( keyPrefix )
      {}

      virtual std::string solvermode () const
      {
        const std::string modes [] = { "dDDI" , "dDFI", "dFFI", "hDDI", "hDFI", "hFFI" };
        int mode = this->parameter_.getEnum(keyPrefix_ + "amgx.mode", modes, 0 );
        return modes[ mode ];
      }

      virtual std::string solverconfig () const
      {
        return this->parameter_.template getValue< std::string >( keyPrefix_ + "amgx.config", "amgxconfig.json");
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

      typedef Fem::PetscLinearOperator< DiscreteFunction, DiscreteFunction > AssembledOperatorType;

      typedef AMGXInverseOperator< DiscreteFunction >  InverseOperatorType;
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

      typedef typename BaseType::SolverDiscreteFunctionType  SolverDiscreteFunctionType;
      typedef typename BaseType::OperatorType                OperatorType;
      typedef typename BaseType::PreconditionerType          PreconditionerType;
      typedef typename BaseType::AssembledOperatorType       AssembledOperatorType;

      /** \brief constructor
       *
       *  \param[in] parameter parameter for the solver
       */
      AMGXInverseOperator ( const SolverParameter& parameter = SolverParameter( Parameter::container() ) )
        : BaseType( parameter )
      {
      }

      void bind( const OperatorType& op )
      {
        BaseType::bind( op );
        const AMGXSolverParameter* param = dynamic_cast< const AMGXSolverParameter* > (&parameter());
        if( param )
        {
          init( *param );
        }
        else
        {
          AMGXSolverParameter newParam( parameter().parameter() );
          init( newParam );
        }
      }

      void unbind()
      {
        amgXSolver_.finalize();
        BaseType :: unbind();
      }

    protected:
      void init( const AMGXSolverParameter& parameter )
      {
        if( assembledOperator_ )
        {
          std::string mode   = parameter.solvermode();
          std::string config = parameter.solverconfig();
          amgXSolver_.initialize(PETSC_COMM_WORLD, mode, config );

          // check that PetscMat was assembled not in block mode
          if( assembledOperator_->blockedMode() )
            DUNE_THROW(InvalidStateException, "AMGXInverseOperator only works with PetscLinearOperator in non-blocked mode!");

          // attach Matrix to linear solver context
          Mat& A = const_cast<Mat &> (assembledOperator_->petscMatrix());

          // set matrix
          amgXSolver_.setA( A );
        }
      }

      int apply( const SolverDiscreteFunctionType& arg, SolverDiscreteFunctionType& dest ) const
      {
        if( !assembledOperator_ )
          DUNE_THROW(NotImplemented,"AMGX solver with matrix free implementations is not supported!");

        // need to have a 'distributed' destination vector for continuous spaces
        if( dest.space().continuous() )
          dest.dofVector().clearGhost();

        // call PETSc solvers, dest = x, arg = rhs
        Vec& x = *dest.petscVec();
        Vec& rhs = *(const_cast< SolverDiscreteFunctionType& > (arg).petscVec());
        amgXSolver_.solve( x, rhs );

        // a continuous solution is 'distributed' so need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

        int iterations = -1;
        // get number of iterations
        amgXSolver_.getIters( iterations );
        return iterations;
      }

    protected:
      mutable AmgXSolver amgXSolver_;
      using BaseType :: assembledOperator_;
      using BaseType :: parameter_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCSOLVER_HH
