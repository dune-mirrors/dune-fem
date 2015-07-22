#ifndef DUNE_FEM_SOLVER_VIENNACL_HH
#define DUNE_FEM_SOLVER_VIENNACL_HH

#define VIENNACL_WITH_EIGEN 
// #define VIENNACL_WITH_OPENCL 
#define VIENNACL_WITH_OPENMP 

#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/vector.hpp>

#include <dune/fem/operator/linear/eigenoperator.hh>

namespace Dune
{

  namespace Fem
  {

    template< class DiscreteFunction >
    class ViennaCLBiCGStabInverseOperator
      : public Fem::Operator< DiscreteFunction, DiscreteFunction >
    {
      typedef Fem::Operator< DiscreteFunction, DiscreteFunction > Base;

    public:
      typedef typename Base::DomainFunctionType DomainFunction;
      typedef typename Base::RangeFunctionType RangeFunction;

      typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;

    private:
      typedef typename OperatorType::MatrixType Matrix;

      typedef typename Matrix::Ttype Field;

    public:
      ViennaCLBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter, bool verbose )
        : matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        viennacl::copy( op.matrix().data(), matrix_ );
      }

      ViennaCLBiCGStabInverseOperator ( const OperatorType &op, double redEps, double absLimit, unsigned int maxIter = std::numeric_limits< unsigned int >::max() )
        : matrix_( op.matrix().rows(), op.matrix().cols() ),
          tag_( absLimit, maxIter )
      {
        viennacl::copy( op.matrix().data(), matrix_ );
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

      typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;

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

      typedef Fem::EigenLinearOperator< RangeFunction, DomainFunction > OperatorType;

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

#endif // #ifndef DUNE_FEM_SOLVER_VIENNACL_HH
