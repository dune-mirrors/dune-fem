#ifndef DUNE_FEM_SOLVER_VIENNACL_HH
#define DUNE_FEM_SOLVER_VIENNACL_HH

#if HAVE_VIENNACL

#if HAVE_EIGEN
#define VIENNACL_WITH_EIGEN
#endif

// OpenCL overrules OpenMP
#if HAVE_OPENCL
//#warning "Using OpneCL"
//#define VIENNACL_WITH_OPENCL
#elif _OPENMP
#error
#define VIENNACL_WITH_OPENMP
#endif

#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/vector.hpp>

#include <dune/fem/solver/inverseoperatorinterface.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>

#include <dune/fem/solver/parameter.hh>

namespace Dune
{

  namespace Fem
  {

    template< class DiscreteFunction, int method = -1 >
    class ViennaCLInverseOperator;

    template< class DiscreteFunction, int method >
    struct ViennaCLInverseOperatorTraits
    {
      typedef DiscreteFunction    DiscreteFunctionType;
      typedef AdaptiveDiscreteFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > SolverDiscreteFunctionType;

      typedef Dune::Fem::Operator< DiscreteFunction, DiscreteFunction > OperatorType;
      typedef OperatorType  PreconditionerType;

      typedef Fem::SparseRowLinearOperator< DiscreteFunction, DiscreteFunction > AssembledOperatorType;

      typedef ViennaCLInverseOperator< DiscreteFunction, method >  InverseOperatorType;
      typedef SolverParameter SolverParameterType;
    };





    template< class DiscreteFunction, int method >
    class ViennaCLInverseOperator
      : public InverseOperatorInterface< ViennaCLInverseOperatorTraits< DiscreteFunction, method > >
    {
      typedef InverseOperatorInterface< ViennaCLInverseOperatorTraits< DiscreteFunction, method > >  BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunction;
      typedef typename BaseType::RangeFunctionType RangeFunction;

      typedef typename DomainFunction :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction  :: DiscreteFunctionSpaceType RangeSpaceType;

#if HAVE_EIGEN
      typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > EigenOperatorType;
#endif

      typedef Fem::SparseRowLinearOperator< RangeFunction, DomainFunction > OperatorType;

      typedef typename BaseType :: SolverDiscreteFunctionType   SolverDiscreteFunctionType;

    private:
      typedef typename OperatorType::MatrixType Matrix;

      typedef typename Matrix::field_type Field;
      typedef typename Matrix::size_type  size_type;

      typedef viennacl::compressed_matrix< Field > ViennaCLMatrix;
      typedef viennacl::vector< Field > ViennaCLVector;

      //typedef viennacl::linalg::ilu0_precond< ViennaCLMatrix > Preconditioner;
      typedef viennacl::linalg::no_precond Preconditioner;
    public:
      ViennaCLInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : ViennaCLInverseOperator( redEps, absLimit, maxIter, verbose )
      {
        bind( op );
      }

      ViennaCLInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : ViennaCLInverseOperator( redEps, absLimit, maxIter, false )
      {
        bind( op );
      }

      ViennaCLInverseOperator ( double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max(),
                                const bool verbose = false,
                                const SolverParameter& parameter = SolverParameter( Parameter::container() ) )
        : BaseType( parameter ),
          absLimit_( absLimit ),
          method_( method < 0 ? parameter.solverMethod({ SolverParameter::gmres, SolverParameter::cg, SolverParameter::bicgstab }) : method )
      {
      }

      ViennaCLInverseOperator ( const SolverParameter& parameter = SolverParameter( Parameter::container() ) )
        : ViennaCLInverseOperator( -1,
            parameter.tolerance(), parameter.maxIterations(),
            parameter.verbose(), parameter )
      {
        assert( parameter.errorMeasure() == 0 );
      }

      void bind( const OperatorType& op )
      {
        BaseType::bind( op );
#if HAVE_EIGEN
        eigenMatrix_ = dynamic_cast<const EigenOperatorType* > (&op);
        if( eigenMatrix_ )
        {
          const auto& matrix = eigenMatrix_->matrix();
          gupMatrix_.resize( matrix.rows(), matrix.cols() );
          viennacl::copy( matrix.data(), gpuMatrix_ );
        }
        else
#endif // HAVE_EIGEN
        if( assembledOperator_ )
        {
          const auto& matrix = assembledOperator_->matrix();
          gpuMatrix_.resize( matrix.rows(), matrix.cols() );
          std::vector< std::map< size_type, Field > > cpuMatrix;
          matrix.fillCSRStorage( cpuMatrix );
          viennacl::copy( cpuMatrix, gpuMatrix_ );
        }
      }

      void unbind()
      {
        BaseType::unbind();
        gpuMatrix_ = ViennaCLMatrix();
#if HAVE_EIGEN
        eigenOperator_ = nullptr;
#endif
      }

      int apply( const SolverDiscreteFunctionType& u, SolverDiscreteFunctionType& w ) const
      {
        viennacl::vector< Field > vclU( u.size() ), vclW( w.size() );
        viennacl::copy( u.dbegin(), u.dend(), vclU.begin() );

        int maxIterations = parameter().maxIterations();
        int iterations = -1;

        //std::cout << "Using ViennaCL " << std::endl;

        if( method_ == SolverParameter::cg )
        {
          Preconditioner ilu0; // ( gpuMatrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::cg_tag tag( absLimit_, maxIterations );
          vclW = viennacl::linalg::solve( gpuMatrix_, vclU, tag, ilu0 );
          iterations = tag.iters();
        }
        else if ( method_ == SolverParameter::bicgstab )
        {
          /*
          // configuration of preconditioner:
          viennacl::linalg::chow_patel_tag chow_patel_ilu_config;
          chow_patel_ilu_config.sweeps(3);       // three nonlinear sweeps
          chow_patel_ilu_config.jacobi_iters(2); // two Jacobi iterations per triangular 'solve' Rx=r
          // create and compute preconditioner:
          typedef viennacl::linalg::chow_patel_ilu_precond< ViennaCLMatrix > Preconditioner;
          //viennacl::linalg::chow_patel_ilu_precond< viennacl::compressed_matrix<ScalarType> > chow_patel_ilu(A, chow_patel_ilu_config);

          typedef viennacl::linalg::block_ilu_precond< ViennaCLMatrix > Preconditioner;
          Preconditioner ilu0( gpuMatrix_, chow_patel_ilu_config );
          */

          //typedef viennacl::linalg::block_ilu_precond< ViennaCLMatrix, viennacl::linalg::ilu0_tag > Preconditioner;
          Preconditioner ilu0; // ( gpuMatrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::bicgstab_tag tag( absLimit_, maxIterations );
          vclW = viennacl::linalg::solve( gpuMatrix_, vclU, tag, ilu0 );
          iterations = tag.iters();
        }
        else if ( method_ == SolverParameter::gmres )
        {
          Preconditioner ilu0; // ( gpuMatrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::gmres_tag tag( absLimit_, maxIterations );
          vclW = viennacl::linalg::solve( gpuMatrix_, vclU, tag, ilu0 );
          iterations = tag.iters();
        }
        else
        {
          DUNE_THROW(NotImplemented,"ViennaCL does not support this solver");
        }

        viennacl::copy( vclW.begin(), vclW.end(), w.dbegin() );
        return iterations ;
      }

    protected:
      ViennaCLMatrix gpuMatrix_;

      using BaseType :: assembledOperator_;
      using BaseType :: iterations_;
      using BaseType :: parameter;
#if HAVE_EIGEN
      const EigenOperatorType* eigenOperator_ = nullptr;
#endif

      double absLimit_;


      int method_;
    };


    //  ViennaCLBiCGStabInverseOperator
    //-----------------------------------

    template< class DiscreteFunction >
    using ViennaCLCGInverseOperator = ViennaCLInverseOperator< DiscreteFunction, SolverParameter :: cg >;


    //  ViennaCLBiCGStabInverseOperator
    //-----------------------------------

    template< class DiscreteFunction >
    using ViennaCLBiCGStabInverseOperator = ViennaCLInverseOperator< DiscreteFunction, SolverParameter :: bicgstab >;


    //  ViennaCLGMResInverseOperator
    //-----------------------------------

    template< class DiscreteFunction >
    using ViennaCLGMResInverseOperator = ViennaCLInverseOperator< DiscreteFunction, SolverParameter :: gmres >;

  } // namespace Fem

} // namespace Dune

#endif

#endif // #ifndef DUNE_FEM_SOLVER_VIENNACL_HH
