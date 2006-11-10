// Emacs should recognise this header as -*- C++ -*-

#ifndef CGHS_BLAS_H
#define CGHS_BLAS_H

// ============================================================================
//
//  CG nach Hestenes und Stiefel
//
//  siehe auch:
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

#include <utility>

// ohne Vorkonditionierer
template< class MATRIX >
inline std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const double *b, double *x, double eps );

template< class MATRIX >
inline std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const double *b, double *x, double eps,
      bool detailed );

template<class CommunicatorType, class MATRIX >
inline int
cghsParallel( CommunicatorType  & comm, unsigned N, const MATRIX &A, const double *b, double *x, double eps,
      bool detailed );


// mit Vorkonditionierer
template< class MATRIX, class PC_MATRIX >
inline std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const PC_MATRIX &C,
      const double *b, double *x, double eps );

template< class MATRIX, class PC_MATRIX >
inline std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const PC_MATRIX &C,
      const double *b, double *x, double eps, bool detailed );


// ============================================================================

#include "cblas.h"
#include "tmpmem.hh"

static OEMTmpMem cghsMem;

template< class CommunicatorType ,class MATRIX >
inline int
cghsParallel( CommunicatorType & comm, unsigned N, const MATRIX &A, const double *b, double *x, double eps,
      bool detailed ) 
{
  if ( N==0 )
    return -1;
#if 1
  cghsMem.resize( 3*N ); 

  double *p = cghsMem.getMem (N);
  double *r = cghsMem.getMem (N);
  double *g = cghsMem.getMem (N);
#else 
  double *g = new double[N];
  double *r = new double[N];
  double *p = new double[N];
#endif
  
  int its=0;
  //double t, tau, sig, rho, gam;
  double t, gam;

  // send and recive buffer for rho,tau and sig
  double * commVal  = new double [3];
  double * commBuff = new double [3];
  
  double & rho_s = commVal[0];
  double & sig_s = commVal[1];
  double & tau_s = commVal[2];
  
  double & rho = commBuff[0];
  double & sig = commBuff[1];
  double & tau = commBuff[2];
  
  double err=eps*eps* comm.globalSum( ddot(N,b,1,b,1) );
  
  mult(A,x,g);
  daxpy(N,-1.,b,1,g,1);
  dscal(N,-1.,g,1);
  dcopy(N,g,1,r,1);
  //comm.communicate();
  double ddo = comm.globalSum( ddot(N,g,1,g,1) );
  while ( ddo>err ) 
  {
    mult(A,r,p);

    rho_s = ddot(N,p,1,p,1) ;
    sig_s = ddot(N,r,1,p,1) ;
    tau_s = ddot(N,g,1,r,1) ;

    comm.globalSumVec ( commVal , 3 , commBuff );
    
    t=tau/sig;
    daxpy(N,t,r,1,x,1);
    daxpy(N,-t,p,1,g,1);
    gam=(t*t*rho-tau)/tau;
    dscal(N,gam,r,1);
    daxpy(N,1.,g,1,r,1);
    if ( detailed )
      std::cout<<"cghs "<<its<<"\t"<<dnrm2(N,g,1)<<std::endl;
    ++its;
    ddo = comm.globalSum( ddot(N,g,1,g,1) );
  }
#if 1
  cghsMem.reset();
#else 
  delete[] g;
  delete[] r;
  delete[] p;
#endif
  delete [] commVal ;
  delete [] commBuff ;
  comm.communicate();
  return its;
}

template< class MATRIX >
inline 
std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const double *b, double *x, double eps,
      bool detailed ) 
{
  if ( N==0 )
  {
    std::cerr << "WARNING: N = 0 in cghs, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }
#if 1
  cghsMem.resize( 3*N ); 

  double *p = cghsMem.getMem (N);
  double *r = cghsMem.getMem (N);
  double *g = cghsMem.getMem (N);
#else 
  double *g = new double[N];
  double *r = new double[N];
  double *p = new double[N];
#endif
  
  int its=0;
  double t, tau, sig, rho, gam;
  double err=eps*eps*ddot(N,b,1,b,1);
  
  mult(A,x,g);
  daxpy(N,-1.,b,1,g,1);
  dscal(N,-1.,g,1);
  dcopy(N,g,1,r,1);
  while ( ddot(N,g,1,g,1)>err ) 
  {
    mult(A,r,p);
    rho=ddot(N,p,1,p,1);
    sig=ddot(N,r,1,p,1);
    tau=ddot(N,g,1,r,1);
    t=tau/sig;
    daxpy(N,t,r,1,x,1);
    daxpy(N,-t,p,1,g,1);
    gam=(t*t*rho-tau)/tau;
    dscal(N,gam,r,1);
    daxpy(N,1.,g,1,r,1);
    if ( detailed )
      std::cout<<"cghs "<<its<<"\t"<<dnrm2(N,g,1)<< std::endl;
    ++its;
  }
  std::pair<int,double> val (its,dnrm2(N,g,1));
  
#if 1
  cghsMem.reset();
#else 
  delete[] g;
  delete[] r;
  delete[] p;
#endif
  return val;
}


template< class MATRIX > 
inline 
std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const double *b, double *x, double eps ) {
  return cghs(N,A,b,x,eps,false);
}


// ============================================================================


template< class MATRIX, class PC_MATRIX >
inline
std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const PC_MATRIX &C,
      const double *b, double *x, double eps, bool detailed ) 
{
  if ( N==0 )
  {
    std::cerr << "WARNING: N = 0 in cghs, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  double *r = new double[N];
  double *d = new double[N];
  double *h = new double[N];

  double *Ad = h;
  int its=0;
  double rh, alpha, beta;
  double err=eps*eps*ddot(N,b,1,b,1);

  mult(A,x,r);
  daxpy(N,-1.,b,1,r,1);
  mult(C,r,d);
  dcopy(N,d,1,h,1);
  rh=ddot(N,r,1,h,1);
  while ( ddot(N,r,1,r,1)>err ) 
  {
    mult(A,d,Ad);
    alpha=rh/ddot(N,d,1,Ad,1);
    daxpy(N,-alpha,d,1,x,1);
    daxpy(N,-alpha,Ad,1,r,1);
    mult(C,r,h);
    beta=1./rh; rh=ddot(N,r,1,h,1); beta*=rh;
    dscal(N,beta,d,1);
    daxpy(N,1.,h,1,d,1);
    if ( detailed )
      cout<<"cghs "<<its<<"\t"<<dnrm2(N,r,1)<<endl;
    ++its;
  }

  std::pair<int,double> val(its,dnrm2(N,r,1));
  
  delete[] r;
  delete[] d;
  delete[] h;

  return val; 
}

template< class MATRIX, class PC_MATRIX > 
inline
std::pair < int , double > 
cghs( unsigned N, const MATRIX &A, const PC_MATRIX &C,
      const double *b, double *x, double eps ) {
  return cghs(N,A,C,b,x,eps,false);
}

// ============================================================================


#endif // CGHS_BLAS_H
