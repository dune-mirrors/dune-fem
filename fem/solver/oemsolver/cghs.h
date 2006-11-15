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
// Modification for parallel use: Robert Kloefkorn
//
// ============================================================================

#include <utility>

// ============================================================================
#include "cblas.h"
#include "tmpmem.hh"

static OEMTmpMem cghsMem;

// CommunicatorType has to fulfill interface of CollectiveCommunication
// see dune-common/common/collectivecommunication.hh
template< class CommunicatorType, class MATRIX >
inline 
std::pair < int , double > 
cghs( const CommunicatorType & comm, 
      unsigned int N, const MATRIX &A, const double *b, double *x, double eps,
      bool detailed ) 
{
  if ( N==0 )
  {
    std::cerr << "WARNING: N = 0 in cghs, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

#ifdef USE_MEMPROVIDER
  cghsMem.resize( 3*N ); 

  double *p = cghsMem.getMem (N);
  double *r = cghsMem.getMem (N);
  double *g = cghsMem.getMem (N);
#else 
  double *g = new double[N];
  double *r = new double[N];
  double *p = new double[N];
#endif
  
  for(size_t k =0 ; k<N; ++k) 
  {
    g[k] = 0.0;
    r[k] = 0.0;
    p[k] = 0.0;
  }
  
  int its=0;
  double t, gam;
  // send and recive buffer for rho,tau and sig
  double commVal[3] = { 0.0,0.0,0.0} ; 
  double * commBuff = (double *) &commVal[0];
  
  double & rho = commVal[0];
  double & sig = commVal[1];
  double & tau = commVal[2];
  
  double bb = ddot(N,b,1,b,1);
  double err=eps*eps* comm.sum( bb );
  
  mult(A,x,g);
  daxpy(N,-1.,b,1,g,1);
  dscal(N,-1.,g,1);
  dcopy(N,g,1,r,1);

  double gg = ddot(N,g,1,g,1);
  double ddo = comm.sum( gg );
  while ( ddo>err ) 
  {
    mult(A,r,p);

    rho=ddot(N,p,1,p,1);
    sig=ddot(N,r,1,p,1);
    tau=ddot(N,g,1,r,1);
    
    comm.sum ( commBuff , 3 );

    t=tau/sig;
    daxpy(N,t,r,1,x,1);
    
    daxpy(N,-t,p,1,g,1);
    gam=(t*t*rho-tau)/tau;
    dscal(N,gam,r,1);
    daxpy(N,1.,g,1,r,1);
    
    gg = ddot(N,g,1,g,1);
    ddo = comm.sum( gg );
    
    if ( detailed && (comm.rank() == 0) )
    {
      std::cout<<"cghs "<<its<<"\t"<<sqrt(ddo)<< std::endl;
    }
    ++its;
  }

  std::pair<int,double> val (its,sqrt(gg));
  
#ifdef USE_MEMPROVIDER
  cghsMem.reset();
#else 
  delete[] g;
  delete[] r;
  delete[] p;
#endif
  
  return val;
}


// ============================================================================

// CommunicatorType has to fulfill interface of CollectiveCommunication
// see dune-common/common/collectivecommunication.hh

template<class CommunicatorType, class MATRIX, class PC_MATRIX >
inline
std::pair < int , double > 
cghs( const CommunicatorType & comm,
      unsigned int N, const MATRIX &A, const PC_MATRIX &C,
      const double *b, double *x, double eps, bool detailed ) 
{
  if ( N==0 )
  {
    std::cerr << "WARNING: N = 0 in cghs, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  // just to remember that parallel is not tested yet 
  assert( comm.size() == 1 );

#ifdef USE_MEMPROVIDER
  cghsMem.resize( 3*N ); 

  double *r = cghsMem.getMem (N);
  double *d = cghsMem.getMem (N);
  double *h = cghsMem.getMem (N);
#else 
  double *r = new double[N];
  double *d = new double[N];
  double *h = new double[N];
#endif
  double *Ad = h;
  int its=0;
  double rh, alpha, beta, rho, sig;
  
  double bb = ddot(N,b,1,b,1);
  double err=eps*eps * comm.sum( bb );

  mult(A,x,r);
  daxpy(N,-1.,b,1,r,1);
  mult(C,r,d);
  dcopy(N,d,1,h,1);
  rh=ddot(N,r,1,h,1);

  double rr = ddot(N,r,1,r,1);
  double ddo = comm.sum( rr );
  while ( ddo>err ) 
  {
    mult(A,d,Ad);
    rho = ddot(N,d,1,Ad,1);
    sig = comm.sum( rho );
    alpha=rh/sig;
    
    daxpy(N,-alpha,d,1,x,1);
    daxpy(N,-alpha,Ad,1,r,1);
    
    mult(C,r,h);
    beta=1./rh; 
    sig = ddot(N,r,1,h,1); 
    rh = comm.sum( sig );

    beta*=rh;
    dscal(N,beta,d,1);
    daxpy(N,1.,h,1,d,1);
    
    rr = ddot(N,r,1,r,1);
    ddo = comm.sum( rr );
    
    if ( detailed && comm.rank() == 0)
    {
        std::cout<<"cghs "<<its<<"\t"<<sqrt(ddo)<< std::endl;
    }
    ++its;
  }

  std::pair<int,double> val(its,sqrt(rr));
  
#ifdef USE_MEMPROVIDER
  cghsMem.reset();
#else 
  delete[] r;
  delete[] d;
  delete[] h;
#endif

  return val; 
}

#endif // CGHS_BLAS_H
