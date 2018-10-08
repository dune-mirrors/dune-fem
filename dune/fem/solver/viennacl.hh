#ifndef DUNE_FEM_SOLVER_VIENNACL_HH
#define DUNE_FEM_SOLVER_VIENNACL_HH

#ifdef HAVE_VIENNACL

#if HAVE_EIGEN
#define VIENNACL_WITH_EIGEN
#endif

#if HAVE_OPENCL
#define VIENNACL_WITH_OPENCL
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
//#include <dune/fem/operator/linear/eigenoperator.hh>

#include <dune/fem/solver/parameter.hh>

namespace Dune
{

  namespace Fem
  {
    template< class DiscreteFunction >
    class ViennaCLInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      typedef typename DomainFunction :: DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunction  :: DiscreteFunctionSpaceType RangeSpaceType;

      //typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;
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
          method_( parameter.krylovMethod() )
      {
      }

      void bind( const OperatorType& op )
      {
        init( op );
      }

      void init( const OperatorType &op )
      {
        matrix_.resize( op.matrix().rows(), op.matrix().rows() );
        std::vector< std::map< size_type, Field > > cpuMatrix;
        op.matrix().fillCSRStorage( cpuMatrix );

        viennacl::copy( cpuMatrix, matrix_ );
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

        if( method_ == SolverParameter::cg )
        {
          Preconditioner ilu0( matrix_, viennacl::linalg::ilu0_tag() );
          viennacl::linalg::cg_tag tag( absLimit_, maxIter_ );
          vclW = viennacl::linalg::solve( matrix_, vclU, tag, ilu0 );
          iterations_ = tag.iters();
        }
        else if ( method_ == SolverParameter::bicgstab )
        {
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

      unsigned int iterations () const { return iterations_; }

    protected:
      ViennaCLMatrix matrix_;

      double absLimit_;
      unsigned int maxIter_;
      mutable unsigned int iterations_;

      int method_;

      mutable std::unique_ptr< ADDomainFunctionType > U_;
      mutable std::unique_ptr< ADRangeFunctionType  > W_;
    };



    template< class DiscreteFunction >
    class ViennaCLBiCGStabInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      //typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;
      typedef Fem::SparseRowLinearOperator< RangeFunction, DomainFunction > OperatorType;

    private:
      typedef typename OperatorType::MatrixType Matrix;

      typedef typename Matrix::field_type Field;
      typedef typename Matrix::size_type  size_type;

    public:
      ViennaCLBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        init( op );
      }

      ViennaCLBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        init( op );
      }

      ViennaCLBiCGStabInverseOperator ( double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max(),
                                        const bool verbose = false )
        : //matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
      }

      void bind( const OperatorType& op )
      {
        init( op );
      }

      void init( const OperatorType &op )
      {
        matrix_.resize( op.matrix().rows(), op.matrix().rows() );
        std::vector< std::map< size_type, Field > > cpuMatrix;
        op.matrix().fillCSRStorage( cpuMatrix );

        viennacl::copy( cpuMatrix, matrix_ );
      }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        viennacl::vector< Field > vclU( u.size() ), vclW( w.size() );
        viennacl::copy( u.dbegin(), u.dend(), vclU.begin() );
        vclW = viennacl::linalg::solve( matrix_, vclU, tag_ );
        viennacl::copy( vclW.begin(), vclW.end(), w.dbegin() );
      }

      unsigned int iterations () const { return tag_.iters(); }

    protected:
      viennacl::compressed_matrix< Field > matrix_;
      viennacl::linalg::bicgstab_tag tag_;
    };



    template< class DiscreteFunction >
    class ViennaCLCGInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      typedef Fem::SparseRowLinearOperator< RangeFunction, DomainFunction > OperatorType;

    private:
      typedef typename OperatorType::MatrixType::MatrixStorageType Matrix;

      typedef typename OperatorType::MatrixType::Ttype Field;

    public:
      ViennaCLCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : op_(op),
          matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        viennacl::copy( op.matrix().data(), matrix_ );
      }

      ViennaCLCGInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : op_(op),
          matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        viennacl::copy( op.matrix().data(), matrix_ );
      }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        std::cout << "starting to solve" << std::endl;
#if 1
        viennacl::vector< Field > vclU( u.size() ), vclW( w.size() );
        viennacl::copy( u.dofVector().array().coefficients(), vclU );
        for (int i=0;i<100;++i)
        {
          std::cout << "starting to solve" << std::endl;
          vclW = viennacl::linalg::solve( matrix_, vclU, tag_ );
        }

        viennacl::copy(  vclW,w.dofVector().array().coefficients() );
#else
        w.dofVector().array().coefficients() = viennacl::linalg::solve(
                 op_.matrix().data(), u.dofVector().array().coefficients(), tag_);
#endif

      }

      unsigned int iterations () const { return tag_.iters(); }

    protected:
      const OperatorType &op_;
      viennacl::compressed_matrix< Field > matrix_;
      viennacl::linalg::cg_tag tag_;
    };



    template< class DiscreteFunction >
    class ViennaCLGMResInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      //typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;

      //typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;
      typedef Fem::SparseRowLinearOperator< RangeFunction, DomainFunction > OperatorType;

    private:
      typedef typename OperatorType::MatrixType Matrix;

      typedef typename Matrix::Ttype Field;

      typedef viennacl::compressed_matrix< Field > ViennaCLMatrix;
      typedef viennacl::vector< Field > ViennaCLVector;

      typedef viennacl::linalg::ilu0_precond< ViennaCLMatrix > Preconditioner;

    public:
      ViennaCLGMResInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        viennacl::copy( op.matrix().data(), matrix_ );
      }

      ViennaCLGMResInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        viennacl::copy( op.matrix().data(), matrix_ );
      }

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        ViennaCLVector vclU( u.size() ), vclW( w.size() );
        Preconditioner ilu0( matrix_, viennacl::linalg::ilu0_tag() );
        viennacl::copy( u.dbegin(), u.dend(), vclU.begin() );
        vclW = viennacl::linalg::solve( matrix_, vclU, tag_, ilu0 );
        viennacl::copy( vclW.begin(), vclW.end(), w.dbegin() );
      }

      unsigned int iterations () const { return tag_.iters(); }

    protected:
      ViennaCLMatrix matrix_;
      viennacl::linalg::gmres_tag tag_;
    };

  } // namespace Fem

} // namespace Dune

#endif

#endif // #ifndef DUNE_FEM_SOLVER_VIENNACL_HH
