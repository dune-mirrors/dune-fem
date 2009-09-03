#ifndef DUNE_UMFPACKSOLVER_HH
#define DUNE_UMFPACKSOLVER_HH

#include <limits>

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#ifdef ENABLE_UMFPACK 
#include <umfpack.h>
#endif

namespace Dune
{

  /** @addtogroup DirectSolver  
      
      In this section implementations of direct solvers 
      for solving linear systems of the from 
      \f$A x = b\f$, where \f$A\f$ is a Mapping or
      Operator and \f$x\f$ and \f$b\f$ are discrete functions 
      (see DiscreteFunctionInterface) can be found. 
   **/


  /** \class UMFPACKOp
   *  \ingroup DirectSolver
   *  \brief UMFPACK direct solver
   */
  template< class DF, class Op >
  struct UMFPACKOp
  : public Operator< typename DF::DomainFieldType, typename DF::RangeFieldType, DF, DF >
  {
    typedef DF DiscreteFunctionType;
    typedef Op OperatorType;

    /** \brief constructor of UMFPACKOp
        \param[in] op Operator to invert 
        \param[in] redEps realative tolerance for residual 
        \param[in] absLimit absolut solving tolerance for residual 
        \param[in] maxIter maximal number of iterations performed 
        \param[in] verbose verbosity 
    */
    UMFPACKOp ( const OperatorType &op, 
                double redEps,
                double absLimit,
                int maxIter,
                bool verbose )
    : op_( op ),
      epsilon_( absLimit ),
      maxIter_( maxIter ),
      verbose_( verbose )
    {}

    UMFPACKOp ( const OperatorType &op,
                double redEps,
                double absLimit,
                int maxIter = std::numeric_limits< int >::max() )
    : op_( op ),
      epsilon_( absLimit ),
      maxIter_( maxIter ),
      verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) )
    {}
                

    void prepare ( const DiscreteFunctionType &, DiscreteFunctionType & ) const
    {}

    void finalize () const
    {}

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    void apply ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) const
    {
      // prepare operator 
      prepare( arg, dest );

#ifdef ENABLE_UMFPACK 
      // call UMF solve method on SparseRowMatrix
      op_.systemMatrix().solveUMF( arg, dest );
#else 
      DUNE_THROW( InvalidStateException, "UMFPACK was not found, reconfigure or use other solver!" );
#endif

      // finalize operator  
      finalize ();
    }

    /** \brief solve the system 
        \param[in] arg right hand side 
        \param[out] dest solution 
    */
    void operator() ( const DiscreteFunctionType &arg, DiscreteFunctionType &dest ) const
    {
      apply( arg, dest );
    }

    void printTexInfo(std::ostream& out) const
    {
      out << "Solver: UMFPACK direct solver ";
      out  << "\\\\ \n";
    }

    double averageCommTime() const 
    {
      return 0.0;
    }

    int iterations() const 
    {
      return 0;
    }

  private:
    // note: the matrix is changed by this operator!
    const OperatorType &op_;
    typename DiscreteFunctionType::RangeFieldType epsilon_;
    int maxIter_;
    bool verbose_ ;
  };

}

#endif // #ifndef DUNE_UMFPACKSOLVER_HH
