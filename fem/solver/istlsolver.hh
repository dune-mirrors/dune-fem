#ifndef DUNE_ISTLSOLVERS_HH 
#define DUNE_ISTLSOLVERS_HH 


#if HAVE_DUNE_ISTL 
//- Dune includes 
#include <dune/fem/operator/common/operator.hh>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

namespace Dune {

  //! Default implementation for the scalar case
  template<class X ,class CommunicatorType>
  class ParaScalarProduct : public ScalarProduct<X>
  {
    const CommunicatorType& comm_;
  public:
    ParaScalarProduct(const CommunicatorType& comm) : comm_(comm) {} 

    //! export types
    typedef X domain_type;
    typedef typename X::field_type field_type;
    
    //! define the category
    enum {category=SolverCategory::sequential};
    
    /*! \brief Dot product of two vectors. 
      It is assumed that the vectors are consistent on the interior+border
      partition.
     */
    virtual field_type dot (const X& x, const X& y)
    {
      field_type val = x*y;
      val = comm_.sum( val );
      return val;
    }

    /*! \brief Norm of a right-hand side vector. 
      The vector must be consistent on the interior+border partition
     */
    virtual double norm (const X& x)
    {
      double val = x.two_norm();
      val = comm_.sum( val ); 
      return val;
    }
  };

  //=====================================================================
  // Implementation for ISTL-matrix based operator
  //=====================================================================

  /*! 
    \brief Adapter to turn a matrix into a linear operator.
    Adapts a matrix to the assembled linear operator interface
  */
  template<class MatrixType, class X, class Y>
  class ParallelMatrixAdapter 
    : public AssembledLinearOperator<MatrixType,X,Y>
  {
  public:
    //! export types
    typedef MatrixType  matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    //! define the category
    enum {category=SolverCategory::sequential};

    //! constructor: just store a reference to a matrix
    ParallelMatrixAdapter (const MatrixType& A) : matrix_(A) {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
      matrix_.mult(x,y);
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
      matrix_.multAdd(alpha,x,y);
    }

    //! get matrix via *
    virtual const MatrixType& getmat () const
    {
      return matrix_;
    }

  private:
    const MatrixType& matrix_;
  };

  
// BICG STAB scheme 
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
      typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;
      typedef typename OperatorImp :: PreconditionMatrixType PreconditionerType; 
      const PreconditionerType& pre = op.preconditionMatrix();
      solve(op.systemMatrix().matrix(),pre,
            arg,dest,arg.space().grid().comm(),reduction,maxIter,verbose);
    }

    template <class MatrixType, 
              class PreconditionerType,
              class DiscreteFunctionImp,
              class CommunicatorType>
    static void solve(const MatrixType & m,
                 const PreconditionerType & preconditioner,
                 const DiscreteFunctionImp & arg,
                 DiscreteFunctionImp & dest,
                 const CommunicatorType& comm,
                 double reduction, int maxIter, bool verbose)
    {
      typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;
      typedef ParallelMatrixAdapter<MatrixType,BlockVectorType,BlockVectorType> MatrixOperatorType;
      MatrixOperatorType mat(const_cast<MatrixType&> (m));

      int verb = (verbose) ? 2 : 0;
        
      ParaScalarProduct<BlockVectorType,CommunicatorType> scp(comm); 
      BiCGSTABSolver<BlockVectorType> solver(mat,scp,
          const_cast<PreconditionerType&> (preconditioner),
          reduction,maxIter,verb);    

      InverseOperatorResult returnInfo;
  
      solver.apply(dest.blockVector(),arg.blockVector(),returnInfo);
    }
  };

public:
  //! ISTL BiCGStab
  ISTLBICGSTABOp(OperatorType & op , double  reduction , double absLimit , 
                int maxIter , bool verbose ) 
    : op_(op), reduction_ ( reduction ) 
    , maxIter_ (maxIter ) , verbose_ ( verbose ) 
  {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    SolverCaller<OperatorType,true>::call(op_,arg,dest,reduction_,maxIter_,verbose_);
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }
}; 


}
#endif

#endif
