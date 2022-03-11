#ifndef DUNE_FEM_FEMPRECONDITIONING_HH
#define DUNE_FEM_FEMPRECONDITIONING_HH

#include <type_traits>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/solver/parameter.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    // FemPreconditioningBase (default non-assembled)
    // --------------------------------------------------
    template< class DFImp, class OperatorImp, int method, bool assembled >
    class FemPreconditioningBase
      : public Operator< DFImp, DFImp >
    {
    public:
      typedef DFImp        DiscreteFunctionType;
      typedef OperatorImp  OperatorType;

      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType;

    public:
      FemPreconditioningBase(const OperatorType &op) {}

      virtual void operator()(const DiscreteFunctionType &u, DiscreteFunctionType &res) const
      {
        DUNE_THROW(NotImplemented,"preconditioning not possible for non-assembled operators");
      }

#if HAVE_DUNE_ISTL
      //! apply for ISTL BlockVectors
      template < class YBlock, class XBlock >
      void applyToISTLBlockVector( const BlockVector< YBlock >& d,
                                   BlockVector< XBlock >& v ) const
      {
        DUNE_THROW(NotImplemented,"preconditioning not possible for non-assembled operators");
      }
#endif

    protected:
      void apply( const DiscreteFunctionType& u, DiscreteFunctionType& res ) const
      {
        DUNE_THROW(NotImplemented,"preconditioning not possible for non-assembled operators");
      }
    };


    // FemPreconditioningBase (assembled version)
    // ----------------------------------------------
    template< class DFImp, class OperatorImp, int method >
    class FemPreconditioningBase< DFImp, OperatorImp, method, true  >
      : public Operator< DFImp, DFImp >
    {
    public:
      typedef DFImp              DiscreteFunctionType;
      typedef OperatorImp        OperatorType;

      typedef typename DiscreteFunctionType :: DofType DofType;
      typedef typename Dune::FieldTraits< DofType >::real_type RealType;
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType;

    protected:
      typedef typename OperatorType::MatrixType MatrixType;

      DiscreteFunctionType diagonalInv_;

      const MatrixType& matrix_;

      const int n_;
      const double w_;

    public:
      FemPreconditioningBase( const OperatorType& assembledOperator,
                              const int n = 1,
                              const double relax = 1.0 )
        : diagonalInv_( "diagInv", assembledOperator.rangeSpace()),
          matrix_( assembledOperator.exportMatrix() ),
          n_( n ),
          w_( (method == SolverParameter::gauss_seidel ) ? 1.0 : relax )
      {
        // extract diagonal elements form matrix object
        assembledOperator.extractDiagonal( diagonalInv_ );

        // make consistent at border dofs
        diagonalInv_.communicate();

        // In general: store 1/diag
        //
        // note: We set near-zero entries to 1 to avoid NaNs. Such entries occur
        //       if DoFs are excluded from matrix setup
        const RealType eps = 16.*std::numeric_limits< RealType >::epsilon();
        const DofIteratorType dend = diagonalInv_.dend();
        for( DofIteratorType dit = diagonalInv_.dbegin(); dit != dend; ++dit )
          *dit = (std::abs( *dit ) < eps ? DofType( 1. ) : DofType( 1. ) / *dit);
      }

      virtual void operator()(const DiscreteFunctionType &u, DiscreteFunctionType &res) const
      {
        apply(u, res);
      }

#if HAVE_DUNE_ISTL
      //! apply for ISTL BlockVectors
      template < class YBlock, class XBlock >
      void applyToISTLBlockVector( const BlockVector< YBlock >& d,
                                   BlockVector< XBlock >& v ) const
      {
        DiscreteFunctionType vTmp("fem-precon::X", diagonalInv_.space(), v );
        DiscreteFunctionType dTmp("fem-precon::Y", diagonalInv_.space(), d );

        // apply stationary iterative preconditioning
        apply( dTmp, vTmp );
      }
#endif

    protected:
      void apply( const DiscreteFunctionType& u, DiscreteFunctionType& v ) const
      {
        v.clear();

        std::unique_ptr< DiscreteFunctionType > xTmp_;
        if constexpr ( method == SolverParameter :: jacobi )
        {
          xTmp_.reset( new DiscreteFunctionType( v ) );
        }

        DiscreteFunctionType& x = ( method == SolverParameter :: jacobi ) ? *xTmp_ : v ;

        for( int i=0; i<n_; ++i )
        {
          matrix_.forwardIterative( diagonalInv_, u, x, v, w_ );

          if constexpr ( method == SolverParameter :: ssor )
          {
            matrix_.backwardIterative( diagonalInv_, u, x, v, w_ );
          }

          // synchronize data
          v.communicate();

          if constexpr ( method == SolverParameter :: jacobi )
          {
            // only needed for the intermediate steps
            if( i < n_-1 )
              x.assign( v );
          }
        }
      }
    };


    // FemPreconditioning
    // ----------------------
    /** \class FemPreconditioning
      *  \ingroup OEMSolver
      *  \brief   Precondtioner, implementing Jacobi, Gauss-Seidel and SOR
      *           works with
      *           - OEM
      *           - ISTL
      *
      *  \param  DFImp     type of the disctete function
      *  \param  Operator  type of the operator (only works for assembled operators)
      */
    template< class DFImp, class Operator, int method>
    class FemPreconditioning
      : public FemPreconditioningBase< DFImp, Operator, method, std::is_base_of< AssembledOperator< DFImp, DFImp >, Operator > :: value >
    {
      typedef FemPreconditioningBase< DFImp, Operator, method, std::is_base_of< AssembledOperator< DFImp, DFImp >, Operator > :: value >
        BaseType;
    public:
      typedef Operator   OperatorType;
      FemPreconditioning(const OperatorType &op, const int n = 1, const double w = 1.0)
        : BaseType( op )
      {}
    };

    template <class DFImp, class Operator>
    using FemJacobiPreconditioning = FemPreconditioning< DFImp, Operator, SolverParameter::jacobi >;

    template <class DFImp, class Operator>
    using FemGaussSeidelPreconditioning = FemPreconditioning< DFImp, Operator, SolverParameter::gauss_seidel >;

    template <class DFImp, class Operator>
    using FemSORPreconditioning = FemPreconditioning< DFImp, Operator, SolverParameter::sor >;

    template <class DFImp, class Operator>
    using FemSSORPreconditioning = FemPreconditioning< DFImp, Operator, SolverParameter::ssor >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DIAGONALPRECONDITIONER_HH
