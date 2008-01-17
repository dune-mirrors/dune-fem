#ifndef DUNE_ISTLSOLVERS_HH 
#define DUNE_ISTLSOLVERS_HH 

#if HAVE_DUNE_ISTL 
//- Dune includes 
#include <dune/fem/operator/common/operator.hh>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/fem/function/common/scalarproducts.hh>

namespace Dune {

  //=====================================================================
  // Implementation for ISTL-matrix based operator
  //=====================================================================

/** @ingroup OEMSolver
    @{
**/
  
/** \brief BICG-stab scheme for block matrices (BCRSMatrix) 
and block vectors (BVector) from dune-istl. */
template <class DiscreteFunctionType, class OperatorType>
class ISTLBICGSTABOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> 
{
private:
  // no const reference, we make const later 
  mutable OperatorType &op_;
  double reduction_;
  int maxIter_;
  bool verbose_ ;

  template <class OperatorImp, bool hasPreconditioning>
  struct SolverCaller
  {
    template <class DiscreteFunctionImp>
    static void call(OperatorImp & op,
                     const DiscreteFunctionImp & arg,
                     DiscreteFunctionImp & dest,
                     double reduction, int maxIter, bool verbose)
    {
      solve(op.systemMatrix(),
            arg,dest,reduction,maxIter,verbose);
    }

    template <class MatrixObjType, 
              class DiscreteFunctionImp>
    static void solve(const MatrixObjType & mObj,
                 const DiscreteFunctionImp & arg,
                 DiscreteFunctionImp & dest,
                 double absLimit, int maxIter, bool verbose)
    {
      typedef typename MatrixObjType :: MatrixAdapterType MatrixAdapterType;
      MatrixAdapterType matrix = mObj.matrixAdapter();
      
      typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

      // verbose only in verbose mode and for rank 0 
      int verb = (verbose && (dest.space().grid().comm().rank() == 0)) ? 2 : 0;
        
      double residuum = matrix.residuum( arg.blockVector(), dest.blockVector());
      double reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

      if( verbose ) 
      {
        std::cout << "ISTL BiCG-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit<< "\n";
      }

      BiCGSTABSolver<BlockVectorType>
        solver(matrix,matrix.scp(),matrix.preconditionAdapter(),
               reduction,maxIter,verb);    

      InverseOperatorResult returnInfo;
  
      solver.apply(dest.blockVector(),arg.blockVector(),returnInfo);
    }
  };

public:
  /** \brief constructor of ISTLBICGSTABOp 
    \param[in] op Mapping describing operator to invert 
    \param[in] redEps reduction epsilon 
    \param[in] absLimit absolut limit of residual (not used here) 
    \param[in] maxIter maximal iteration steps 
    \param[in] verbose verbosity 

    \note ISTL BiCG-stab only uses the relative reduction.
  */
  ISTLBICGSTABOp(OperatorType & op , double  reduction , double absLimit , 
                int maxIter , bool verbose ) 
    : op_(op), reduction_ ( absLimit ) 
    , maxIter_ (maxIter ) , verbose_ ( verbose ) 
  {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }


  /** \brief solve the system 
      \param[in] arg right hand side 
      \param[out] dest solution 
  */
  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    SolverCaller<OperatorType,true>::call(op_,arg,dest,reduction_,maxIter_,verbose_);
  }


  /** \brief solve the system 
      \param[in] arg right hand side 
      \param[out] dest solution 
  */
  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }
}; 


//////////////////////////////////////////////////////////////////
//
//  ISTL CG Solver 
//
//////////////////////////////////////////////////////////////////
/** \brief BICG-stab scheme for block matrices (BCRSMatrix) 
and block vectors (BVector) from dune-istl. */
template <class DiscreteFunctionType, class OperatorType>
class ISTLCGOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> 
{
private:
  // no const reference, we make const later 
  mutable OperatorType &op_;
  const double absLimit_;
  int maxIter_;
  bool verbose_ ;

  template <class OperatorImp, bool hasPreconditioning>
  struct SolverCaller
  {
    template <class DiscreteFunctionImp>
    static void call(OperatorImp & op,
                     const DiscreteFunctionImp & arg,
                     DiscreteFunctionImp & dest,
                     double absLimit, int maxIter, bool verbose)
    {
      solve(op.systemMatrix(),
            arg,dest,absLimit,maxIter,verbose);
    }

    template <class MatrixObjType, 
              class DiscreteFunctionImp>
    static void solve(const MatrixObjType & mObj,
                 const DiscreteFunctionImp & arg,
                 DiscreteFunctionImp & dest,
                 double absLimit, int maxIter, bool verbose)
    {
      typedef typename MatrixObjType :: MatrixAdapterType MatrixAdapterType;
      MatrixAdapterType matrix = mObj.matrixAdapter();
      
      typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;

      // verbose only in verbose mode and for rank 0 
      int verb = (verbose && (dest.space().grid().comm().rank() == 0)) ? 2 : 0;
        
      double residuum = matrix.residuum( arg.blockVector(), dest.blockVector());
      double reduction = (residuum > 0) ? absLimit/ residuum : 1e-3;

      if( verbose ) 
      {
        std::cout << "ISTL CG-Solver: reduction: " << reduction << ", residuum: " << residuum << ", absolut limit: " << absLimit<< "\n";
      }

      CGSolver<BlockVectorType> 
        solver(matrix,matrix.scp(),matrix.preconditionAdapter(),
          reduction,maxIter,verb);    

      InverseOperatorResult returnInfo;
  
      solver.apply(dest.blockVector(),arg.blockVector(),returnInfo);
    }
  };

public:
  /** \brief constructor of ISTLBICGSTABOp 
    \param[in] op Mapping describing operator to invert 
    \param[in] redEps reduction epsilon 
    \param[in] absLimit absolut limit of residual (not used here) 
    \param[in] maxIter maximal iteration steps 
    \param[in] verbose verbosity 

    \note ISTL BiCG-stab only uses the relative reduction.
  */
  ISTLCGOp(OperatorType & op , double  reduction , double absLimit , 
           int maxIter , bool verbose ) 
    : op_(op), absLimit_ ( absLimit ) 
    , maxIter_ (maxIter ) , verbose_ ( verbose ) 
  {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }


  /** \brief solve the system 
      \param[in] arg right hand side 
      \param[out] dest solution 
  */
  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    SolverCaller<OperatorType,true>::call(op_,arg,dest,absLimit_,maxIter_,verbose_);
  }


  /** \brief solve the system 
      \param[in] arg right hand side 
      \param[out] dest solution 
  */
  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }
}; 
///@}

} // end namespace Dune 
#endif // end HAVE_DUNE_ISTL

#endif
