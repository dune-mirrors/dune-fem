#ifndef DUNE_ISTLSOLVERS_HH 
#define DUNE_ISTLSOLVERS_HH 


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

  template<class X, class Y>
  class EmptyPreconditioner : public Preconditioner<X,Y> {
  public:

    EmptyPreconditioner () {}
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    enum {
      //! \brief The category the precondtioner is part of.
      category=SolverCategory::sequential};

    /*! \brief Prepare the preconditioner. 

    A solver solves a linear operator equation A(x)=b by applying 
        one or several steps of the preconditioner. The method pre()
        is called before before the first apply operation. 
        x and b are right hand side and solution vector of the linear
        system. It may. e.g., scale the system, allocate memory or
        compute a (I)LU decomposition.
      Note: The ILU decomposition could also be computed in the constructor
        or with a separate method of the derived method if several
        linear systems with the same matrix are to be solved.

        \param x The left hand side of the equation.
        \param b The right hand side of the equation.
    */
    virtual void pre (X& x, Y& b) {
      //x = b; 
    }

    /*! \brief Apply one step of the preconditioner to the system A(v)=d. 

        On entry v=0 and d=b-A(x) (although this might not be 
        computed in that way. On exit v contains the update, i.e
        one step computes \f$ v = M^{-1} d \f$ where \f$ M \f$ is the
        approximate inverse of the operator \f$ A \f$ characterizing 
        the preconditioner.
        \param[out] v The update to be computed
        \param d The current defect.
    */
    virtual void apply (X& v, const Y& d)
    {
      v = d;
    }

    /*! \brief Clean up.

    This method is called after the last apply call for the
        linear system to be solved. Memory may be deallocated safely
        here. x is the solution of the linear equation.

        \param x The right hand side of the equation.
    */
    virtual void post (X& x) {
    }

    // every abstract base class has a virtual destructor
    virtual ~EmptyPreconditioner () {}
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
      matrix_.usmv(alpha,x,y);
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
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  template <class OperatorImp, bool hasPreconditioning>
  struct SolverCaller
  {
    template <class DiscreteFunctionImp>
    static void call(OperatorImp & op,
                     const DiscreteFunctionImp & arg,
                     DiscreteFunctionImp & dest,
                     double eps, int maxIter, bool verbose)
    {
      typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;
      typedef typename OperatorImp :: PreconditionMatrixType PreconMatrix;
      
      if( op.hasPreconditionMatrix() ) 
      {
        solve(op.systemMatrix().matrix(),op.preconditionMatrix(),
              arg,dest,arg.space().grid().comm(),eps,maxIter,verbose);
      }
      else 
      {
        EmptyPreconditioner<BlockVectorType,BlockVectorType> preconditioner;
        solve(op.systemMatrix().matrix(),preconditioner,
            arg,dest,arg.space().grid().comm(),eps,maxIter,verbose);
      }
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
                 double eps, int maxIter, bool verbose)
    {
      typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;
      typedef ParallelMatrixAdapter<MatrixType,BlockVectorType,BlockVectorType> MatrixOperatorType;
      MatrixOperatorType mat(const_cast<MatrixType&> (m));

      int verb = (verbose) ? 2 : 0;
        
      ParaScalarProduct<BlockVectorType,CommunicatorType> scp(comm); 
      BiCGSTABSolver<BlockVectorType> solver(mat,scp,
          const_cast<PreconditionerType&> (preconditioner),
          eps,maxIter,verb);    

      InverseOperatorResult returnInfo;
  
      solver.apply(dest.blockVector(),arg.blockVector(),returnInfo);
    }
  };

public:
  ISTLBICGSTABOp(OperatorType & op , double  redEps , double absLimit , 
                int maxIter , bool verbose ) 
    : op_(op), epsilon_ ( sqrt(absLimit) ) 
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
    SolverCaller<OperatorType,true>::call(op_,arg,dest,epsilon_,maxIter_,verbose_);
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }
}; 


}
#endif
