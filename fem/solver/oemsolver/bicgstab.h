// -*- C++ -*-

#ifndef BICGSTAB_BLAS_H
#define BICGSTAB_BLAS_H

// ============================================================================
//
//  BICGstab
//
//  siehe
//  Ashby, Manteuffel, Saylor
//     A taxononmy for conjugate gradient methods
//     SIAM J Numer Anal 27, 1542-1568 (1990)
//
//  oder:
//  Willy D"orfler:
//     Orthogonale Fehlermethoden
//
//                                                 ----------------------------
//                                                 Christian Badura, Mai 1998
//   
//  Modifications and parallelization: Robert Kloefkorn  
//
// ============================================================================

#include <utility>

// ============================================================================

#include <iostream>
#include <assert.h>
#include "cblas.h"

// ============================================================================
#include "tmpmem.hh"

static OEMTmpMem bicgMem;

template<bool usePC, 
         class CommunicatorType, 
         class MATRIX ,
         class PC_MATRIX > 
inline
std::pair<int,double> 
bicgstab_algo( const CommunicatorType & comm,
    unsigned int N, const MATRIX &A, const PC_MATRIX & C,
	  const double *rhs, double *x, double eps, bool detailed ) 
{
  if(N == 0) 
  {
    std::cerr << "WARNING: N = 0 in bicgstab, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  typedef Mult<MATRIX,PC_MATRIX,usePC> MultType; 

  double * tmp = 0;
#ifdef USE_MEMPROVIDER
  int memSizeFactor = 7;
  if( usePC ) ++memSizeFactor;
  bicgMem.resize ( memSizeFactor * N ); 
  
  double *r   = bicgMem.getMem( N );
  double *t   = bicgMem.getMem( N );
  double *Ad  = bicgMem.getMem( N );
  double *u   = bicgMem.getMem( N );
  double *s   = bicgMem.getMem( N );
  double *d   = bicgMem.getMem( N );
  double *rT  = bicgMem.getMem( N );
  if( usePC ) tmp = bicgMem.getMem( N );
#else 
  double *rT  = new double[N];
  double *d   = new double[N];
  double *s   = new double[N];
  double *u   = new double[N];
  double *Ad  = new double[N];
  double *t   = new double[N];
  double *r   = new double[N];
  if( usePC ) tmp = new double[N];
#endif

  const double bicgeps = 1e-40;

  // f"ur's Abbruchkriterium (*)  -- r enth"alt immer das Residuum r=Ax-b
  unsigned int its=0;

  double st, tt, rTr, rTh; 
  
  double rTAd, alpha, beta, omega;

  double err=eps*eps;
  double bb = 0.0;

  bool rightPreCon = MultType :: mult_pc(A,C,x,r,tmp);

  // if pc matrix, recalc rhs 
  if( usePC && (!rightPreCon) )
  {
    mult(C,rhs,tmp);
    daxpy(N,-1.,tmp,1,r,1);
    bb = MultType::ddot(A,tmp,tmp); 
  }
  else 
  {
    daxpy(N,-1.,rhs,1,r,1);
    bb = MultType::ddot(A,rhs,rhs); 
  }

  err *= bb ; 
  
  dcopy(N,r,1,d,1);
  dcopy(N,d,1,s,1);
  dcopy(N,s,1,rT,1);
 
  rTr = MultType::ddot(A,r,r);
  rTh = MultType::ddot(A,rT,s);
 
  while( rTr>err ) 
  {
    // do multiply 
    MultType :: mult_pc(A,C,d,Ad,tmp);
    rTAd = MultType::ddot(A,rT,Ad);
   
    // if no error, return 
    if( fabs(rTAd) <= bicgeps ) break;

    alpha = rTh/rTAd;

    daxpy(N,-alpha,Ad,1,r,1);
    daxpy(N,-alpha,Ad,1,s,1);
    
    // do multiply 
    MultType :: mult_pc(A,C,s,t,tmp);

    daxpy(N,1.,t,1,u,1);
    dscal(N,alpha,u,1);

    st = MultType::ddot(A,s,t);
    tt = MultType::ddot(A,t,t);

    if ( fabs(st)<bicgeps || fabs(tt)<bicgeps )
    {
      omega = 0.;
    }
    else
    {
      omega = st/tt;
    }

    daxpy(N,-omega,t,1,r,1);
    daxpy(N,-alpha,d,1,x,1);
    daxpy(N,-omega,s,1,x,1);
    
    daxpy(N,-omega,t,1,s,1);
    beta=(alpha/omega)/rTh; 
    
    rTh = MultType::ddot(A,rT,s); 
    rTr = MultType::ddot(A,r,r);

    beta*=rTh;
    dscal(N,beta,d,1);
    daxpy(N,1.,s,1,d,1);
    daxpy(N,-beta*omega,Ad,1,d,1);

    if ( detailed && (comm.rank() == 0) )
    {
      std::cout<<"bicgstab "<<its<<"\t  tol: " << err << "   err: "<<rTr<< std::endl;
    }
    ++its;
  }

  // if right preconditioning then do back solve 
  MultType :: back_solve(N,C,x,tmp);
       
  std::pair<int,double> val (its,sqrt(rTr));

#ifdef USE_MEMPROVIDER
  bicgMem.reset();
#else
  delete[] r;
  delete[] rT;
  delete[] d;
  delete[] s;
  delete[] u;
  delete[] Ad;
  delete[] t;
  if( usePC ) delete [] tmp;
#endif
  
  return val;
}


// bicgstab with pc matrix 
template<class CommunicatorType,
         class MATRIX > 
inline
std::pair<int,double> 
bicgstab( const CommunicatorType & comm,
    unsigned int N, const MATRIX &A,
	  const double *b, double *x, double eps, bool verbose ) 
{
  return bicgstab_algo<false>(comm,N,A,A,b,x,eps,verbose);
}

// bicgstab with pc matrix 
template<class CommunicatorType,
         class MATRIX,
         class PC_MATRIX> 
inline
std::pair<int,double> 
bicgstab( const CommunicatorType & comm,
    unsigned int N, const MATRIX &A, const PC_MATRIX & C,
	  const double *b,double *x, double eps, bool verbose ) 
{
  return bicgstab_algo<true>(comm,N,A,C,b,x,eps,verbose);
}

#endif // BICGSTAB_BLAS_H
