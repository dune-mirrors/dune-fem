#ifndef __DUNE_OEMSOLVER_HH__
#define __DUNE_OEMSOLVER_HH__

#include "../../operator/common/operator.hh"

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
  
#include "bicgstab.h"
#include "cghs.h"
#include "gmres.h"
#include "bicgsq.h"
  
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
    typedef typename DiscreteFunctionType::FunctionSpaceType FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeFieldType Field;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    typedef typename DiscreteFunctionType::ConstDofIteratorType ConstDofIteratorType;
    typedef typename FunctionSpaceType::GridType GridType;

    // prepare operator 
    prepare ( arg, dest );

    int size = arg.getFunctionSpace().size();
    
    OEMSolver::cghs
      (size,op_.systemMatrix(),arg.leakPointer(),dest.leakPointer(),epsilon_,verbose_);

    // finalize operator  
    finalize ();
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};


// cummunicate the unkown of the CG Method 
template <class Writer , class Reader>
class CommunicateVector
{
  Writer & dataWriter_;
  Reader & dataReader_;
public:
  CommunicateVector ( Writer & w, Reader & r ) : dataWriter_(w), dataReader_(r) {}

  template <class ObjectStreamType, class EntityType>
  void inlineData ( ObjectStreamType & str, EntityType & en )
  {
    assert(false);
  }

  template <class ObjectStreamType, class EntityType>
  void xtractData ( ObjectStreamType & str, EntityType & en )
  {
    assert(false);
  }

  template <class ObjectStreamType, class EntityType>
  void scatter ( ObjectStreamType & str, EntityType & en )
  {
    std::pair < ObjectStreamType * , EntityType * > p (&str,&en);
    dataWriter_.apply( p );
  }

  template <class ObjectStreamType, class EntityType>
  void gather ( ObjectStreamType & str, EntityType & en )
  {
    std::pair < ObjectStreamType * , EntityType * > p (&str,&en);
    dataReader_.apply( p );
  }
};

// LocalCommType = CommunicateVector
template <class GridCommType , class LocalCommType>
class CommunicateCG
{
  GridCommType & gc_;
  LocalCommType & comm_;
public:
  CommunicateCG ( GridCommType & gc , LocalCommType & comm )
    : gc_(gc), comm_(comm) {}

  void communicate ()
  {
    gc_.communicate(comm_);
  }

  double globalSum( double val )
  {
    return gc_.globalSum(val);
  }

  void globalSumVec ( double * send, int s , double * recv)
  {
    gc_.globalSum(send,s,recv);
  }
};

//! Parallel version of above CG scheme
template <class CommunicatorType , class DiscreteFunctionType, class OperatorType>
class OEMCGOpParallel : public Operator<
      typename DiscreteFunctionType::DomainFieldType,
      typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType,DiscreteFunctionType> {

private:
  // no const reference, we make const later 
  OperatorType &op_;
  typename DiscreteFunctionType::RangeFieldType epsilon_;
  int maxIter_;
  bool verbose_ ;
  CommunicatorType & comm_;
public:

  OEMCGOpParallel( CommunicatorType & comm, OperatorType & op , double  redEps , double absLimit , int maxIter , bool verbose ) :
        comm_(comm), op_(op), epsilon_ ( absLimit ) ,
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
    typedef typename DiscreteFunctionType::FunctionSpace FunctionSpaceType;
    typedef typename FunctionSpaceType::RangeField Field;
    typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;
    typedef typename DiscreteFunctionType::ConstDofIteratorType ConstDofIteratorType;
    typedef typename FunctionSpaceType::GridType GridType;

    // prepare operator 
    prepare ( arg, dest );

    int size = arg.getFunctionSpace().size();
    
    OEMSolver::cghsParallel
      (comm_,size,op_.systemMatrix(),arg.leakPointer(),dest.leakPointer(),epsilon_,verbose_);

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

    int size = arg.getFunctionSpace().size();

    
    //int iter = 
    OEMSolver::bicgstab
      (size,op_.systemMatrix(),arg.leakPointer(),dest.leakPointer(),epsilon_,verbose_);

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

    int size = arg.getFunctionSpace().size();

    
    //int iter = 
    OEMSolver::bicgsq
      (size,op_.systemMatrix(),arg.leakPointer(),dest.leakPointer(),epsilon_,verbose_);

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

    int size = arg.getFunctionSpace().size();
    int inner = (size > 1000) ? 1000 : size;

    OEMSolver::gmres
      (inner,size,op_.systemMatrix(),arg.leakPointer(),dest.leakPointer(),epsilon_,verbose_);

    // finalize operator  
    finalize ();
  }

  void operator ()( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
  {
    apply(arg,dest);
  }

};
 
} // end namespace Dune 

#endif
