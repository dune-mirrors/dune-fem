#ifndef DUNE_OEMSOLVER_HH
#define DUNE_OEMSOLVER_HH

//- system includes 
#include <utility>

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/fem/operator/common/operator.hh>

//- local includes 
#include "preconditioning.hh"

#include "../odesolver.hh"

namespace OEMSolver 
{

//! this method is called from all solvers and is only a wrapper
//! this method is mainly from SparseRowMatrix 
template <class MatrixImp, class VectorType>
void mult(const MatrixImp & m, const VectorType * x, VectorType * ret)
{
  // call multOEM of the matrix 
  m.multOEM(x,ret);
}

//! mult method when given pre conditioning matrix 
template <class Matrix , class PC_Matrix , bool >
struct Mult
{
  typedef void mult_t(const Matrix &A,
                      const PC_Matrix & C, 
                      const double *arg,
                      double *dest , 
                      double * tmp);
  
  static void mult_pc (const Matrix &A, const PC_Matrix & C, const double *arg,
        double *dest , double * tmp)
  {
    assert( tmp );

    // call mult of Matrix A 
    mult(A,arg,tmp);
    // call mult of Matrix PC
    mult(C,tmp,dest);
  }
};

//! mult method when no pre conditioning matrix 
template <class Matrix>
struct Mult<Matrix,Matrix,false>
{
  typedef void mult_t(const Matrix &A,
                      const Matrix & C, 
                      const double *arg,
                      double *dest , 
                      double * tmp);
  
  static void mult_pc(const Matrix &A, const Matrix & C, const double *arg ,
                      double *dest , double * tmp)
  {
    // tmp has to be 0
    assert( tmp == 0 );
    // C is just a fake 
    assert( &A == &C );
    // call mult of Matrix A 
    mult(A,arg,dest);
  }
};

#define USE_MEMPROVIDER   
#include "bicgstab.h"
#include "cghs.h"
#include "gmres.h"
#include "bicgsq.h"
#undef USE_MEMPROVIDER
  
} // end namespace OEMSolver 


namespace Dune 
{

// CG scheme after Hestenes and Stiefel
template <class DiscreteFunctionType, class OperatorType>
class OEMCGOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> {

private:
  // no const reference, we make const later 
  OperatorType &op_;
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  typedef std::pair < int , double > ReturnValueType;

  template <class OperatorImp, bool hasPreconditioning> 
  struct SolverCaller 
  {
    template <class DiscreteFunctionImp> 
    static ReturnValueType call(OperatorImp & op, 
                     const DiscreteFunctionImp & arg, 
                     DiscreteFunctionImp & dest, 
                     double eps, bool verbose)
    {
      // use communication class of grid
      // see dune-common/common/collectivecommunication.hh 
      // for interface 
      int size = arg.space().size();

      if(op.hasPreconditionMatrix())
      {
        return OEMSolver::cghs(arg.space().grid().comm(),
                   size,op.systemMatrix(),op.preconditionMatrix(),
                   arg.leakPointer(),dest.leakPointer(),eps,verbose);
      }
      else 
      {
        return OEMSolver::cghs(arg.space().grid().comm(),
                  size,op.systemMatrix(),
                  arg.leakPointer(),dest.leakPointer(),eps,verbose);
      }
    }
  };

  //! without any preconditioning 
  template <class OperatorImp> 
  struct SolverCaller<OperatorImp,false> 
  {
    template <class DiscreteFunctionImp> 
    static ReturnValueType call(OperatorImp & op, 
                     const DiscreteFunctionImp & arg, 
                     DiscreteFunctionImp & dest, 
                     double eps, bool verbose)
    {
      // use communication class of grid
      // see dune-common/common/collectivecommunication.hh 
      // for interface 
      int size = arg.space().size();
      return OEMSolver::cghs(arg.space().grid().comm(),
                size,op.systemMatrix(),
                arg.leakPointer(),dest.leakPointer(),eps,verbose);
    }
  };

public:

  OEMCGOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose ) :
        op_(op), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_ ( verbose ) {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    // prepare operator 
    prepare ( arg, dest );

    ReturnValueType val = 
      SolverCaller<OperatorType,
                   // check wheter operator has precondition methods 
                   // to enable preconditioning derive your operator from 
                   // OEMSolver::PreconditionInterface
                   Conversion<OperatorType, OEMSolver::PreconditionInterface > ::exists >::
                     // call solver, see above 
                     call(op_,arg,dest,epsilon_,verbose_);

    if(arg.space().grid().comm().rank() == 0)
    {
      std::cout << "OEM-CG: " << val.first << " iterations! Error: " << val.second << "\n";
    }

    // finalize operator  
    finalize ();
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }
};


// BICG STAB scheme 
template <class DiscreteFunctionType, class OperatorType>
class OEMBICGSTABOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> {

private:
  // no const reference, we make const later 
  OperatorType &op_;
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  typedef std::pair < int , double > ReturnValueType;
  
  template <class OperatorImp, bool hasPreconditioning> 
  struct SolverCaller 
  {
    template <class DiscreteFunctionImp> 
    static ReturnValueType call(OperatorImp & op, 
                     const DiscreteFunctionImp & arg, 
                     DiscreteFunctionImp & dest, 
                     double eps, bool verbose)
    {
      int size = arg.space().size();
      if(op.hasPreconditionMatrix())
      {
        return OEMSolver::bicgstab(arg.space().grid().comm(),
                  size,op.systemMatrix(),op.preconditionMatrix(),
                  arg.leakPointer(),dest.leakPointer(),eps,verbose);
      }
      else 
      {
        return OEMSolver::bicgstab(arg.space().grid().comm(),
                  size,op.systemMatrix(),
                  arg.leakPointer(),dest.leakPointer(),eps,verbose);
      }
    }
  };

  //! without any preconditioning 
  template <class OperatorImp> 
  struct SolverCaller<OperatorImp,false> 
  {
    template <class DiscreteFunctionImp> 
    static ReturnValueType call(OperatorImp & op, 
                     const DiscreteFunctionImp & arg, 
                     DiscreteFunctionImp & dest, 
                     double eps, bool verbose)
    {
      int size = arg.space().size();
      return OEMSolver::bicgstab(arg.space().grid().comm(),
                size,op.systemMatrix(),
                arg.leakPointer(),dest.leakPointer(),eps,verbose);
    }
  };

public:

  OEMBICGSTABOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose ) :
        op_(op), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_ ( verbose ) {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

    // prepare operator 
    prepare ( arg, dest );

    ReturnValueType val = 
      SolverCaller<OperatorType,
                   // check wheter operator has precondition methods 
                   // to enable preconditioning derive your operator from 
                   // OEMSolver::PreconditionInterface
                   Conversion<OperatorType, OEMSolver::PreconditionInterface > ::exists >::
                     // call solver, see above 
                     call(op_,arg,dest,epsilon_,verbose_);
    
    if(arg.space().grid().comm().rank() == 0)
    {
      std::cout << "OEM-BICGstab: " << val.first << " iterations! Error: " << val.second << "\n";
    }

    // finalize operator  
    finalize ();
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};

////////////////////////////////
// BICG SQ scheme 
////////////////////////////////
template <class DiscreteFunctionType, class OperatorType>
class OEMBICGSQOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> {

private:
  // no const reference, we make const later 
  OperatorType &op_;
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

public:

  OEMBICGSQOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose ) :
        op_(op), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_ ( verbose ) {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

    // prepare operator 
    prepare ( arg, dest );

    int size = arg.space().size();

    int iter = OEMSolver::bicgsq(size,op_.systemMatrix(),
        arg.leakPointer(),dest.leakPointer(),epsilon_,verbose_);

    std::cout << "OEM-BICGGsq: " << iter << " iterations!\n";
    // finalize operator  
    finalize ();
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};


template <class DiscreteFunctionType, class OperatorType>
class OEMGMRESOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> {

private:
  // no const reference, we make const later 
  OperatorType &op_;
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  typedef std::pair < int , double > ReturnValueType;
  
  template <class OperatorImp, bool hasPreconditioning> 
  struct SolverCaller 
  {
    template <class DiscreteFunctionImp> 
    static ReturnValueType call(OperatorImp & op, 
                     const DiscreteFunctionImp & arg, 
                     DiscreteFunctionImp & dest, 
                     int inner, double eps, bool verbose)
    {
      int size = arg.space().size();
      if(op.hasPreconditionMatrix())
      {
        return OEMSolver::gmres(arg.space().grid().comm(),
                 inner,size,op.systemMatrix(),op.preconditionMatrix(),
                 arg.leakPointer(),dest.leakPointer(),eps,verbose);
      }
      else 
      {
        return OEMSolver::gmres(arg.space().grid().comm(),
                 inner,size,op.systemMatrix(),
                 arg.leakPointer(),dest.leakPointer(),eps,verbose);
      }
    }
  };

  // without any preconditioning 
  template <class OperatorImp> 
  struct SolverCaller<OperatorImp,false>
  {
    template <class DiscreteFunctionImp> 
    static ReturnValueType call(OperatorImp & op, 
                     const DiscreteFunctionImp & arg, 
                     DiscreteFunctionImp & dest, 
                     int inner, double eps, bool verbose)
    {
      int size = arg.space().size();
      return OEMSolver::gmres(arg.space().grid().comm(),
               inner,size,op.systemMatrix(),
               arg.leakPointer(),dest.leakPointer(),eps,verbose);
    }
  };

public:

  OEMGMRESOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose ) :
        op_(op), epsilon_ ( absLimit ) ,
        maxIter_ (maxIter ) , verbose_ ( verbose ) {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;

    // prepare operator 
    prepare ( arg, dest );

    int size = arg.space().size();
    int inner = (size > 20) ? 20 : size;

    ReturnValueType val = 
      SolverCaller<OperatorType,
                   // check wheter operator has precondition methods 
                   // to enable preconditioning derive your operator from 
                   // OEMSolver::PreconditionInterface
                   Conversion<OperatorType, OEMSolver::PreconditionInterface > ::exists >::
                     // call solver, see above 
                     call(op_,arg,dest,inner,epsilon_,verbose_);

    if(arg.space().grid().comm().rank() == 0)
    {
      std::cout << "OEM-GMRES: " << val.first << " iterations! Error: " << val.second << "\n";
    }

    // finalize operator  
    finalize ();
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};

/////////////////////////////////////////////////////////////////
//
//  GMRES Version of Dennis code
//
/////////////////////////////////////////////////////////////////
//! \brief GMRES implementation from Dennis D.
template <class DiscreteFunctionType, class OperatorType>
class GMRESOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> 
{

private:
  // solver 
  mutable DuneODE::GMRES solver_;
  // wrapper to fit interface of GMRES operator 
  mutable DuneODE::SolverInterfaceImpl<OperatorType> op_; 
  
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  typedef std::pair < int , double > ReturnValueType;
  
public:
  GMRESOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose )
      : solver_(DuneODE::Communicator::instance(),20)
      , op_(op), epsilon_ ( absLimit ) 
      , maxIter_ (maxIter ) , verbose_ ( verbose ) 
  {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  //! solve the system 
  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    // prepare operator 
    prepare ( arg, dest );

    int size = arg.space().size();
    op_.setSize(size);

    solver_.set_tolerance(epsilon_);
    solver_.set_max_number_of_iterations(size);

    if(verbose_)
    {
      solver_.IterativeSolver::set_output(std::cout);
      solver_.DynamicalObject::set_output(std::cout);
    }

    // note argument and destination are toggled 
    solver_.solve(op_, dest.leakPointer() , arg.leakPointer() );

    // finalize operator  
    finalize ();
  }

  //! solve the system 
  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};
 
template <class DiscreteFunctionType, class OperatorType>
class FGMRESOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> 
{

private:
  template <class SolverType, bool hasPreconditioning> 
  struct SolverCaller 
  {
    template <class OperatorImp, class PreConMatrix, class DiscreteFunctionImp> 
    static void  call(SolverType & solver, 
               OperatorImp & op, 
               const PreConMatrix & pm, 
               const DiscreteFunctionImp & arg, 
               DiscreteFunctionImp & dest) 
    {
      DuneODE::SolverInterfaceImpl<PreConMatrix> pre(pm); 
      solver.set_preconditioner(pre);

      // note argument and destination are toggled 
      solver.solve(op, dest.leakPointer() , arg.leakPointer() );
    }
  };

  // without any preconditioning 
  template <class SolverType> 
  struct SolverCaller<SolverType,false>
  {
    template <class OperatorImp, class PreConMatrix, class DiscreteFunctionImp> 
    static void  call(SolverType & solver, 
               OperatorImp & op, 
               const PreConMatrix & pm, 
               const DiscreteFunctionImp & arg, 
               DiscreteFunctionImp & dest)
    {
      // note argument and destination are toggled 
      solver.solve(op, dest.leakPointer() , arg.leakPointer() );
    }
  };

  // solver 
  mutable DuneODE::FGMRES solver_;
  // wrapper to fit interface of FGMRES operator 
  typedef DuneODE::FGMRES SolverType;
  mutable DuneODE::SolverInterfaceImpl<OperatorType> op_; 
  mutable OperatorType & orgOp_;
  
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  typedef std::pair < int , double > ReturnValueType;
  
public:
  FGMRESOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose )
      : solver_(DuneODE::Communicator::instance(),20)
      , op_(op), orgOp_(op) , epsilon_ ( absLimit ) 
      , maxIter_ (maxIter ) , verbose_ ( verbose ) 
  {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  //! solve the system 
  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    // prepare operator 
    prepare ( arg, dest );

    int size = arg.space().size();
    op_.setSize(size);

    solver_.set_tolerance(epsilon_);
    solver_.set_max_number_of_iterations(size);

    if(verbose_)
    {
      solver_.IterativeSolver::set_output(std::cout);
      solver_.DynamicalObject::set_output(std::cout);
    }

    if( orgOp_.hasPreconditionMatrix())
    {
      SolverCaller<SolverType,
                   // check wheter operator has precondition methods 
                   // to enable preconditioning derive your operator from 
                   // OEMSolver::PreconditionInterface
                   Conversion<OperatorType, OEMSolver::PreconditionInterface > ::exists >::
                     // call solver, see above 
                     call(solver_,op_,orgOp_.preconditionMatrix(), arg,dest);
      
    }
    else 
    {
      assert(false);
    }

    // finalize operator  
    finalize ();
  }

  //! solve the system 
  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};
 
/////////////////////////////////////////////////////////////////
//
//  BICGstab Version of Dennis code
//
/////////////////////////////////////////////////////////////////
//! \brief GMRES implementation from Dennis D.
template <class DiscreteFunctionType, class OperatorType>
class BICGSTABOp : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> 
{

private:
  // solver 
  mutable DuneODE::BICGSTAB solver_;
  // wrapper to fit interface of GMRES operator 
  mutable DuneODE::SolverInterfaceImpl<OperatorType> op_; 
  
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;

  typedef std::pair < int , double > ReturnValueType;
  
public:
  BICGSTABOp( OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose )
      : solver_(DuneODE::Communicator::instance())
      , op_(op), epsilon_ ( absLimit ) 
      , maxIter_ (maxIter ) , verbose_ ( verbose ) 
  {
  }

  void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
  {
  }

  void finalize () const
  {
  }

  //! solve the system 
  void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    // prepare operator 
    prepare ( arg, dest );

    int size = arg.space().size();
    op_.setSize(size);

    solver_.set_tolerance(epsilon_);
    solver_.set_max_number_of_iterations(size);

    if(verbose_)
    {
      solver_.IterativeSolver::set_output(std::cout);
      solver_.DynamicalObject::set_output(std::cout);
    }

    // note argument and destination are toggled 
    solver_.solve(op_, dest.leakPointer() , arg.leakPointer() );

    // finalize operator  
    finalize ();
  }

  //! solve the system 
  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};
 
} // end namespace Dune 

#endif
