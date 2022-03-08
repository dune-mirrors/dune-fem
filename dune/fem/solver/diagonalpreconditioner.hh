#ifndef DUNE_FEM_DIAGONALPRECONDITIONER_HH
#define DUNE_FEM_DIAGONALPRECONDITIONER_HH

#include <type_traits>

#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/operator/common/operator.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bvector.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    // DiagonalPreconditionerBase (default non-assembled)
    // --------------------------------------------------
    template< class DFImp, class OperatorImp, bool assembled >
    class DiagonalPreconditionerBase
      : public Operator< DFImp, DFImp >
    {
    public:
      typedef DFImp        DiscreteFunctionType;
      typedef OperatorImp  OperatorType;

      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType;

    public:
      DiagonalPreconditionerBase(const OperatorType &op) {}

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


    // DiagonalPreconditionerBase (assembled version)
    // ----------------------------------------------
    template< class DFImp, class OperatorImp >
    class DiagonalPreconditionerBase< DFImp, OperatorImp, true  >
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

    public:
      DiagonalPreconditionerBase( const OperatorType& assembledOperator )
        : diagonalInv_( "diag-preconditioning", assembledOperator.rangeSpace()),
          matrix_( assembledOperator.exportMatrix() )

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
        DiscreteFunctionType vTmp("diag-precon::X", diagonalInv_.space(), v );
        DiscreteFunctionType dTmp("diag-precon::Y", diagonalInv_.space(), d );

        // apply 1/diagonal
        apply( dTmp, vTmp );
      }
#endif


    protected:
      void apply( const DiscreteFunctionType& u, DiscreteFunctionType& res ) const
      {
        /*
        matrix_.sor( diagonalInv_, u, res, 1.0 );
        */

        ConstDofIteratorType uIt     = u.dbegin();
        ConstDofIteratorType diagInv = diagonalInv_.dbegin();

        const DofIteratorType resEnd = res.dend();

        // apply 1/diagonal
        for(DofIteratorType resIt = res.dbegin();
            resIt != resEnd; ++ resIt, ++diagInv, ++ uIt )
        {
          assert( diagInv != diagonalInv_.dend() );
          assert( uIt     != u.dend() );
          (*resIt) = (*uIt) * (*diagInv);
        }
      }

    };


    // DiagonalPreconditioner
    // ----------------------
    /** \class DiagonalPreconditioner
      *  \ingroup OEMSolver
      *  \brief   Precondtioner, multiplies with inverse of the diagonal
      *           works with
      *           - OEM
      *           - ISTL (choose jacobi, iteration = 1)
      *
      *  \param  DFImp     type of the disctete function
      *  \param  Operator  type of the operator (only works for assembled operators)
      */
    template< class DFImp, class Operator>
    class DiagonalPreconditioner
      : public DiagonalPreconditionerBase< DFImp, Operator, std::is_base_of< AssembledOperator< DFImp, DFImp >, Operator > :: value >
    {
      typedef DiagonalPreconditionerBase< DFImp, Operator, std::is_base_of< AssembledOperator< DFImp, DFImp >, Operator > :: value >
        BaseType;
    public:
      typedef Operator   OperatorType;
      DiagonalPreconditioner(const OperatorType &op)
        : BaseType( op )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DIAGONALPRECONDITIONER_HH
