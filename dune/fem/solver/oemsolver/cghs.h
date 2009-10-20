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
template< bool usePC,
          class CommunicatorType, 
          class MATRIX , 
          class PC_MATRIX >
inline 
std::pair < int , double > 
cghs_algo( const CommunicatorType & comm, 
      unsigned int N, const MATRIX &A, const PC_MATRIX& C,
      const double *b, double *x, double eps,
      bool detailed ) 
{
  if ( N==0 )
  {
    std::cerr << "WARNING: N = 0 in cghs, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  // define type of multiplication 
  typedef Mult<MATRIX,PC_MATRIX,usePC> MultType;

  double* tmp = 0;

#ifdef USE_MEMPROVIDER
  cghsMem.resize( (usePC) ? 4*N : 3*N ); 

  double *p = cghsMem.getMem (N);
  double *r = cghsMem.getMem (N);
  double *g = cghsMem.getMem (N);
  if( usePC ) tmp = cghsMem.getMem(N);
#else 
  double *g = new double[N];
  double *r = new double[N];
  double *p = new double[N];
  if( usePC ) tmp = new double [N];
#endif
  
  for(size_t k =0 ; k<N; ++k) 
  {
    g[k] = 0.0;
    r[k] = 0.0;
    p[k] = 0.0;
  }
  
  int its=0;
  double t, gam;
  double rho, sig, tau;
  
  double err= eps * eps * MultType :: ddot(A,b,b);
  
  // apply first multiplication 
  MultType :: mult_pc(A,C,x,g,tmp);
    
  daxpy(N,-1.,b,1,g,1);
  dscal(N,-1.,g,1);
  dcopy(N,g,1,r,1);

  double gg = MultType :: ddot(A,g,g);

  while ( gg > err ) 
  {
    // apply multiplication 
    MultType :: mult_pc(A,C,r,p,tmp);

    rho = MultType :: ddot(A,p,p);
    sig = MultType :: ddot(A,r,p);
    tau = MultType :: ddot(A,g,r);
    
    t=tau/sig;
    daxpy(N,t,r,1,x,1);
    
    daxpy(N,-t,p,1,g,1);
    gam=(t*t*rho-tau)/tau;
    dscal(N,gam,r,1);
    daxpy(N,1.,g,1,r,1);
    
    gg = MultType :: ddot(A,g,g);
    
    if ( detailed && (comm.rank() == 0) )
    {
      std::cout<<"cghs "<<its<<"\t"<<sqrt(gg)<< std::endl;
    }
    ++its;
  }

  // if right preconditioning then do back solve 
  MultType :: back_solve(N,C,x,tmp);

  std::pair<int,double> val (its,sqrt(gg));
  
#ifdef USE_MEMPROVIDER
  cghsMem.reset();
#else 
  delete[] g;
  delete[] r;
  delete[] p;
  if( tmp ) delete [] tmp;
#endif
  
  return val;
}

// cghs with preconditioning 
template<class CommunicatorType, class MATRIX, class PC_MATRIX >
inline
std::pair < int , double > 
cghs( const CommunicatorType & comm,
      unsigned int N, const MATRIX &A, const PC_MATRIX &C,
      const double *b, double *x, double eps, bool detailed ) 
{
  return cghs_algo<true> (comm,N,A,C,b,x,eps,detailed);
}

// cghs without preconditioning 
template<class CommunicatorType, class MATRIX >
inline
std::pair < int , double > 
cghs( const CommunicatorType & comm,
      unsigned int N, const MATRIX &A,
      const double *b, double *x, double eps, bool detailed ) 
{
  return cghs_algo<false> (comm,N,A,A,b,x,eps,detailed);
}

#endif // CGHS_BLAS_H
