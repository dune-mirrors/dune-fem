#ifndef DUNE_FEM_SOLVER_VIENNACL_HH
#define DUNE_FEM_SOLVER_VIENNACL_HH

#if HAVE_VIENNACL

#if HAVE_EIGEN
#define VIENNACL_WITH_EIGEN
#endif

#if HAVE_OPENCL
//#define VIENNACL_WITH_OPENCL
#endif

#if _OPENMP
#define VIENNACL_WITH_OPENMP
#endif

#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/vector.hpp>

#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>

#include <dune/fem/solver/parameter.hh>

namespace Dune
{

  namespace Fem
  {
    template< class DiscreteFunction, int method = -1 >
    class ViennaCLInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      typedef typename DomainFunction :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction  :: DiscreteFunctionSpaceType RangeSpaceType;

#if HAVE_EIGEN
      typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > EigenOperatorType;
#endif

      typedef Fem::SparseRowLinearOperator< RangeFunction, DomainFunction > OperatorType;

      typedef AdaptiveDiscreteFunction< DomainSpaceType > ADDomainFunctionType;
      typedef AdaptiveDiscreteFunction< DomainSpaceType > ADRangeFunctionType;

    private:
      typedef typename OperatorType::MatrixType Matrix;

      typedef typename Matrix::field_type Field;
      typedef typename Matrix::size_type  size_type;

      typedef viennacl::compressed_matrix< Field > ViennaCLMatrix;
      typedef viennacl::vector< Field > ViennaCLVector;

      typedef viennacl::linalg::ilu0_precond< ViennaCLMatrix > Preconditioner;
    public:
      ViennaCLInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : ViennaCLInverseOperator( redEps, absLimit, maxIter, verbose )
      {
        init( op );
      }

      ViennaCLInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : ViennaCLInverseOperator( redEps, absLimit, maxIter, false )
      {
        init( op );
      }

      ViennaCLInverseOperator ( double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max(),
                                const bool verbose = false,
                                const SolverParameter& parameter = SolverParameter( Parameter::container() ) )
        : absLimit_( parameter.linReductionParameter() ),
          maxIter_( maxIter ),
          iterations_( 0 ),
          method_( method < 0 ? parameter.krylovMethod() : method )
      {
      }

      void bind( const OperatorType& op )
      {
        init( op );
      }

    protected:
      void init( const OperatorType &op )
      {
        matrix_.resize( op.matrix().rows(), op.matrix().rows() );
        std::vector< std::map< size_type, Field > > cpuMatrix;
        op.matrix().fillCSRStorage( cpuMatrix );

        viennacl::copy( cpuMatrix, matrix_ );
      }

    public:
#if HAVE_EIGEN
      void bind( const EigenOperatorType& op )
      {
        matrix_.resize( op.matrix().rows(), op.matrix().cols() );
        viennacl::copy( op.matrix().data(), matrix_ );
      }
#endif

      void unbind()
      {
        matrix_ = ViennaCLMatrix();
        U_.reset();
        W_.reset();
      }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        apply( u, w );
      }

      template <class DImpl, class RImpl>
      void operator() ( const DiscreteFunctionInterface< DImpl >& u, DiscreteFunctionInterface< RImpl > &w ) const
      {
        apply( u, w );
      }

      void apply( const ADDomainFunctionType& u, ADRangeFunctionType& w ) const
      {
        viennacl::vector< Field > vclU( u.size() ), vclW( w.size() );
        viennacl::copy( u.dbegin(), u.dend(), vclU.begin() );

        //std::cout << "Using ViennaCL " << std::endl;

        if( method_ == SolverParameter::cg )
        {
          Preconditioner ilu0( matrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::cg_tag tag( absLimit_, maxIter_ );
          vclW = viennacl::linalg::solve( matrix_, vclU, tag, ilu0 );
          iterations_ = tag.iters();
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
          Preconditioner ilu0( matrix_, chow_patel_ilu_config );
          */

          //typedef viennacl::linalg::block_ilu_precond< ViennaCLMatrix, viennacl::linalg::ilu0_tag > Preconditioner;
          Preconditioner ilu0( matrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::bicgstab_tag tag( absLimit_, maxIter_ );
          vclW = viennacl::linalg::solve( matrix_, vclU, tag, ilu0 );
          iterations_ = tag.iters();
        }
        else if ( method_ == SolverParameter::gmres )
        {
          Preconditioner ilu0( matrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::gmres_tag tag( absLimit_, maxIter_ );
          vclW = viennacl::linalg::solve( matrix_, vclU, tag, ilu0 );
          iterations_ = tag.iters();
        }
        else
        {
          DUNE_THROW(NotImplemented,"ViennaCL does not support this solver");
        }

        viennacl::copy( vclW.begin(), vclW.end(), w.dbegin() );
      }

      template <class DImpl, class RImpl>
      void apply( const DiscreteFunctionInterface< DImpl >& u, DiscreteFunctionInterface< RImpl > &w ) const
      {
        if( ! U_ )
        {
          U_.reset( new ADDomainFunctionType("viennacl::u", u.space() ) );
        }

        if( ! W_ )
        {
          W_.reset( new ADRangeFunctionType("viennacl::w", w.space() ) );
        }

        U_->assign( u );

        apply( *U_, *W_ );

        w.assign( *W_ );
      }

      int iterations () const { return iterations_; }

    protected:
      ViennaCLMatrix matrix_;

      double absLimit_;
      unsigned int maxIter_;
      mutable int iterations_;

      int method_;

      mutable std::unique_ptr< ADDomainFunctionType > U_;
      mutable std::unique_ptr< ADRangeFunctionType  > W_;
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
