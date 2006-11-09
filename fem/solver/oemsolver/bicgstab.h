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
// ============================================================================



template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps );


template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps, bool detailed );



template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX &C,
	  const double *b, double *x, double eps );


template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX &C,
	  const double *b, double *x, double eps, bool detailed );


// ============================================================================

#include <iostream>
#include <assert.h>
#include "cblas.h"

// ============================================================================
#include "tmpmem.hh"

static OEMTmpMem bicgMem;

template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps, bool detailed ) 
{

#if 1
  bicgMem.resize ( 7*N ); 
  
  double *r   = bicgMem.getMem( N );
  double *t   = bicgMem.getMem( N );
  double *Ad  = bicgMem.getMem( N );
  double *u   = bicgMem.getMem( N );
  double *h   = bicgMem.getMem( N );
  double *d   = bicgMem.getMem( N );
  double *rT  = bicgMem.getMem( N );
#else 
  double *rT  = new double[N];
  double *d   = new double[N];
  double *h   = new double[N];
  double *u   = new double[N];
  double *Ad  = new double[N];
  double *t   = new double[N];
  double *r   = new double[N];
#endif
  
  // f"ur's Abbruchkriterium (*)  -- r enth"alt immer das Residuum r=Ax-b
 
  double *s   = h;
  double rTh, rTAd, rTr, alpha, beta, omega, st, tt;
  unsigned int its=0;

  double err=eps*eps*ddot(N,b,1,b,1);

  mult(A,x,r);
  daxpy(N,-1.,b,1,r,1);
  dcopy(N,r,1,d,1);
  dcopy(N,d,1,h,1);
  dcopy(N,h,1,rT,1);
  
  //std::cout << "Start defekt = " << ddot(N,rT,1,rT,1) << "\n";
  assert( ddot(N,rT,1,rT,1)>1e-40 );
  rTh=ddot(N,rT,1,h,1);
  rTr=ddot(N,r,1,r,1);
  while ( rTr>err ) {
    mult(A,d,Ad);
    rTAd=ddot(N,rT,1,Ad,1);
    assert( fabs(rTAd)>1e-40 );
    alpha=rTh/rTAd;
    daxpy(N,-alpha,Ad,1,r,1);
    dcopy(N,h,1,s,1);
    daxpy(N,-alpha,Ad,1,s,1);
    mult(A,s,t);
    daxpy(N,1.,t,1,u,1);
    dscal(N,alpha,u,1);
    st=ddot(N,s,1,t,1);
    tt=ddot(N,t,1,t,1);
    if ( fabs(st)<1e-40 || fabs(tt)<1e-40 )
      omega = 0.;
    else
      omega = st/tt;
    daxpy(N,-omega,t,1,r,1);
    daxpy(N,-alpha,d,1,x,1);
    daxpy(N,-omega,s,1,x,1);
    dcopy(N,s,1,h,1);
    daxpy(N,-omega,t,1,h,1);
    beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
    dscal(N,beta,d,1);
    daxpy(N,1.,h,1,d,1);
    daxpy(N,-beta*omega,Ad,1,d,1);
    rTr=ddot(N,r,1,r,1);
    if ( detailed )
      std::cout<<"bicgstab "<<its<<"\t"<<sqrt(rTr)<<std::endl;
    ++its;
  }

#if 1
  bicgMem.reset();
#else
  delete[] r;
  delete[] rT;
  delete[] d;
  delete[] h;
  delete[] u;
  delete[] Ad;
  delete[] t;
#endif
  
  return its;
}



template< class MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A,
	  const double *b, double *x, double eps ) {
  return bicgstab(N,A,b,x,eps,false);
}

// ============================================================================

template< class MATRIX, class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX & C,
	  const double *rhs, double *x, double eps, bool detailed ) 
{

#if 1
  bicgMem.resize ( 9*N ); 
  
  double *r   = bicgMem.getMem( N );
  double *t   = bicgMem.getMem( N );
  double *Ad  = bicgMem.getMem( N );
  double *u   = bicgMem.getMem( N );
  double *h   = bicgMem.getMem( N );
  double *d   = bicgMem.getMem( N );
  double *rT  = bicgMem.getMem( N );
  double *newRhs = bicgMem.getMem( N );
  double *tmp = bicgMem.getMem( N );
#else 
  double *rT  = new double[N];
  double *d   = new double[N];
  double *h   = new double[N];
  double *u   = new double[N];
  double *Ad  = new double[N];
  double *t   = new double[N];
  double *r   = new double[N];
  double *newRhs = new double[N];
  double * tmp = new double[N];
#endif

  for(register int i=0; i<N; ++i) newRhs[i] = 0.0;
  // create new rhs 
  mult(C,rhs,newRhs);
  const double * b = newRhs;
  
  // f"ur's Abbruchkriterium (*)  -- r enth"alt immer das Residuum r=Ax-b
 
  double *s   = h;
  double rTh, rTAd, rTr, alpha, beta, omega, st, tt;
  unsigned int its=0;

  double err=eps*eps*ddot(N,b,1,b,1);

  // mult(A,x,r);
  mult_pc(A,C,x,r,tmp);
  daxpy(N,-1.,b,1,r,1);
  dcopy(N,r,1,d,1);
  dcopy(N,d,1,h,1);
  dcopy(N,h,1,rT,1);
  
  //std::cout << "Start defekt = " << ddot(N,rT,1,rT,1) << "\n";
  assert( ddot(N,rT,1,rT,1)>1e-40 );
  rTh=ddot(N,rT,1,h,1);
  rTr=ddot(N,r,1,r,1);
  while ( rTr>err ) {
    //mult(A,d,Ad);
    mult_pc(A,C,d,Ad,tmp);
    rTAd=ddot(N,rT,1,Ad,1);
    assert( fabs(rTAd)>1e-40 );
    alpha=rTh/rTAd;
    daxpy(N,-alpha,Ad,1,r,1);
    dcopy(N,h,1,s,1);
    daxpy(N,-alpha,Ad,1,s,1);
    //mult(A,s,t);
    mult_pc(A,C,s,t,tmp);
    daxpy(N,1.,t,1,u,1);
    dscal(N,alpha,u,1);
    st=ddot(N,s,1,t,1);
    tt=ddot(N,t,1,t,1);
    if ( fabs(st)<1e-40 || fabs(tt)<1e-40 )
      omega = 0.;
    else
      omega = st/tt;
    daxpy(N,-omega,t,1,r,1);
    daxpy(N,-alpha,d,1,x,1);
    daxpy(N,-omega,s,1,x,1);
    dcopy(N,s,1,h,1);
    daxpy(N,-omega,t,1,h,1);
    beta=(alpha/omega)/rTh; rTh=ddot(N,rT,1,h,1); beta*=rTh;
    dscal(N,beta,d,1);
    daxpy(N,1.,h,1,d,1);
    daxpy(N,-beta*omega,Ad,1,d,1);
    rTr=ddot(N,r,1,r,1);
    if ( detailed )
      std::cout<<"bicgstab "<<its<<"\t"<<sqrt(rTr)<<std::endl;
    ++its;
  }

#if 1
  bicgMem.reset();
#else
  delete[] r;
  delete[] rT;
  delete[] d;
  delete[] h;
  delete[] u;
  delete[] Ad;
  delete[] t;
  delete[] newRhs;
  delete[] tmp;
#endif
  
  return its;
}



template< class MATRIX , class PC_MATRIX > inline
int
bicgstab( unsigned N, const MATRIX &A, const PC_MATRIX & C,
	  const double *b, double *x, double eps ) {
  return bicgstab(N,A,C,b,x,eps,false);
}



// ============================================================================

#endif // BICGSTAB_BLAS_H
