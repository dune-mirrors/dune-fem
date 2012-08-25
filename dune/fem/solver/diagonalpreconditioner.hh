#ifndef DUNE_FEM_DIAGONALPRECONDITIONER_HH
#define DUNE_FEM_DIAGONALPRECONDITIONER_HH

#include <dune/common/static_assert.hh>

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
    template< class DFImp, class AssembledOperator >
    class DiagonalPreconditionerBase< DFImp, AssembledOperator, true >
      : public Operator< DFImp, DFImp >
    {
    public:   
      typedef DFImp              DiscreteFunctionType;
      typedef AssembledOperator  OperatorType;
      
      typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;  
      typedef typename DiscreteFunctionType :: ConstDofIteratorType ConstDofIteratorType;  
        
    protected:
      DiscreteFunctionType diagonalInv_;

    public:
      DiagonalPreconditionerBase( const OperatorType& assembledOperator )
        : diagonalInv_( "diag-preconditioning", assembledOperator.systemMatrix().domainSpace() )
      {
        // estract diagonal elements form matrix object 
        assembledOperator.systemMatrix().extractDiagonal( diagonalInv_ );

        // make consistent at border dofs  
        diagonalInv_.communicate();

        // get dof type 
        typedef typename DiscreteFunctionType :: DofType DofType;

        // In general: store 1/diag 
        //
        // note: set zero entries to 1, this can happen when dofs 
        // are excluded from the matrix setup
        const DofIteratorType diagEnd = diagonalInv_.dend();
        for( DofIteratorType diagInv = diagonalInv_.dbegin();
             diagInv != diagEnd; ++ diagInv ) 
        {
          // get dof entry 
          const DofType& dofii = (*diagInv);
          DofType dof;
          // if dof is zero, store 1 to avoid NaN 
          if( std::abs( dofii ) < 1e-14 )
          {
            dof = 1;
          }
          else 
          {
            dof = 1.0 / dofii;
          }
        }
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
      : public DiagonalPreconditionerBase< DFImp, Operator, Operator::assembled > 
    {
      typedef DiagonalPreconditionerBase< DFImp, Operator, Operator::assembled > BaseType;
    public:
      typedef Operator   OperatorType;
      DiagonalPreconditioner(const OperatorType &op) 
        : BaseType( op )
      {}
    };

  } // namespace Fem 

#if DUNE_FEM_COMPATIBILITY  
  // put this in next version 1.4 

  using Fem :: DiagonalPreconditioner ;
#endif // DUNE_FEM_COMPATIBILITY

} // namespace Dune 

#endif // #ifndef DUNE_FEM_DIAGONALPRECONDITIONER_HH
