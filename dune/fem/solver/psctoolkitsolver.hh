#ifndef DUNE_FEM_PSCTOOLKITSOLVER_HH
#define DUNE_FEM_PSCTOOLKITSOLVER_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/fem/solver/inverseoperatorinterface.hh>
#include <dune/fem/solver/parameter.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/function/adaptivediscretefunction.hh>

#if HAVE_PSCTOOLKIT
// PSCTOOLKIT C bindings for krylov solvers.
#include <psc_krylov_cbind.h>
#endif

namespace Dune
{

  namespace Fem
  {

    //=====================================================================
    // Implementation of PSCTOOLKIT solver wrapper using SparseRowMatrix
    //=====================================================================

    /** @ingroup OEMSolver
        @{
    **/

    struct PSCToolkitSolverParameter : public LocalParameter< SolverParameter, PSCToolkitSolverParameter >
    {
      typedef LocalParameter< SolverParameter, PSCToolkitSolverParameter >  BaseType;

    public:
      using BaseType :: parameter;
      using BaseType :: keyPrefix;

      PSCToolkitSolverParameter( const ParameterReader& parameter = Parameter::container() )
        : BaseType( parameter )
      {}


      PSCToolkitSolverParameter( const std::string &keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : BaseType( keyPrefix, parameter )
      {}

      PSCToolkitSolverParameter( const SolverParameter& sp )
        : PSCToolkitSolverParameter( sp.keyPrefix(), sp.parameter() )
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



    // PSCTOOLKITSolver
    // --------------

    template< class DiscreteFunction >
    class PSCTOOLKITInverseOperator;

    template< class DiscreteFunction >
    struct PSCTOOLKITInverseOperatorTraits
    {
      typedef DiscreteFunction    DiscreteFunctionType;
      typedef PetscDiscreteFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > SolverDiscreteFunctionType;

      typedef Dune::Fem::Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef OperatorType  PreconditionerType;

      typedef Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction > AssembledOperatorType;

      typedef PSCToolkitInverseOperator< DiscreteFunction >  InverseOperatorType;

      typedef PSCToolkitSolverParameter SolverParameterType;
    };




    /** \brief PSCToolkit solver context */
    template< class DF >
    class PSCToolkitInverseOperator
      : public InverseOperatorInterface< PSCToolkitInverseOperatorTraits< DF > >
    {
      typedef PSCToolkitInverseOperatorTraits< DF >  Traits;
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
      PSCToolkitInverseOperator ( const PSCToolkitSolverParameter& parameter = PSCToolkitSolverParameter() )
        : BaseType( parameter )
      {
      }

      PSCToolkitInverseOperator ( const OperatorType &op,
                                  const PSCToolkitSolverParameter & parameter = PSCToolkitSolverParameter() )
        : PSCToolkitInverseOperator ( parameter )
      {
        bind( op );
      }

      PSCToolkitInverseOperator ( const OperatorType &op, PreconditionerType &preconditioner,
                                  const PSCToolkitSolverParameter & parameter = PSCToolkitSolverParameter() )
        : PSCToolkitInverseOperator( parameter )
      {
        bind( op, preconditioner );
      }

      PSCToolkitInverseOperator ( const PSCToolkitInverseOperator& other )
        : PSCToolkitInverseOperator( other.parameter() )
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
#if HAVE_PSCTOOLKIT
        // cleanup
#endif
        BaseType :: unbind();
      }

    protected:
      void init( const PSCToolkitSolverParameter& parameter )
      {
        if( assembledOperator_ )
        {
          std::string mode   = parameter.solvermode();
          std::string config = parameter.solverconfig();
#if HAVE_PSCTOOLKIT
          {
            // TODO copy matrix and call solver

            // attach Matrix to linear solver context
            // Mat& A = const_cast<Mat &> (assembledOperator_->exportMatrix());
          }
#else
          DUNE_THROW(InvalidStateException,"PSCToolkit solver not found during cmake config. Please reconfigure!");
#endif
        }
      }

      int apply( const SolverDiscreteFunctionType& arg, SolverDiscreteFunctionType& dest ) const
      {
        if( !assembledOperator_ )
          DUNE_THROW(NotImplemented,"PSCToolkit solver with matrix free implementations is not supported!");

        int iterations = -1;
#if HAVE_PSCTOOLKIT
#else
        DUNE_THROW(InvalidStateException,"PSCToolkit solver not found during cmake config. Please reconfigure!");
#endif
        return iterations;
      }

    protected:
#if HAVE_PSCTOOLKIT
#endif
      using BaseType :: assembledOperator_;
      using BaseType :: parameter_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PSCTOOLKITSOLVER_HH
