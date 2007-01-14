#ifndef DUNE_ISTLSOLVERS_HH 
#define DUNE_ISTLSOLVERS_HH 

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

namespace Dune {

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
      solve(op.systemMatrix().matrix(),arg,dest,eps,maxIter,verbose);
    }

    template <class MatrixType, class DiscreteFunctionImp>
    static void solve(const MatrixType & m,
                 const DiscreteFunctionImp & arg,
                 DiscreteFunctionImp & dest,
                 double eps, int maxIter, bool verbose)
    {
      //int size = arg.space().size();
      //if(op.hasPreconditionMatrix())
      {
        //int numIterations = 10;
        //double relaxationFactor = 1.1;
        typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;
        //SeqILUn<MatrixType,BlockVectorType,BlockVectorType> 
        //  preconditioner(m,numIterations,relaxationFactor); 
        
        EmptyPreconditioner<BlockVectorType,BlockVectorType> preconditioner;
        typedef MatrixAdapter<MatrixType,BlockVectorType,BlockVectorType>
          MatrixOperatorType;
        MatrixOperatorType mat(const_cast<MatrixType&> (m));

        BiCGSTABSolver<BlockVectorType> solver(mat,preconditioner,
                                               eps,maxIter,2);    

        InverseOperatorResult returnInfo;
        solver.apply(dest.blockVector(),arg.blockVector(),returnInfo);
      }
    }
  };



public:
  ISTLBICGSTABOp(OperatorType & op , double  redEps , double absLimit , 
                int maxIter , bool verbose ) :
        op_(op), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_ ( verbose ) 
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
