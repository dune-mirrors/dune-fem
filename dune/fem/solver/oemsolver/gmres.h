#ifndef GMRES_BLAS_H
#define GMRES_BLAS_H

// ============================================================================
//
//  GMRES nach Saad, Schultz
//     GMRES: a generalized minimal residual algorithm for solving nonsymmetric
//     linear systems
//     SIAM J Sci Stat Comput 7, 856-869 (1986)
//
//                                                 ----------------------------
//                                                 Christian Badura, Mai 1998
//
// ============================================================================

#include <utility>
#include "cblas.h"
#include "tmpmem.hh"

static OEMTmpMem gmresMem;


template<bool usePC ,
         class CommunicatorType, 
         class Matrix , 
         class PC_Matrix >
inline
std::pair<int,double> 
gmres_algo (const CommunicatorType & comm,
       int m, int n, const Matrix &A, const PC_Matrix & C, 
       const double *b , double *x, double eps,
       bool detailed ) 
{
  if ( n<=0 )
  {
    std::cerr << "WARNING: n = " << n << " in gmres_pc, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  // to make sure it's not used for parallel, because not working yet
  if( comm.size() > 1 ) 
  {
    std::cerr << "gmres not working for parallel runs yet! file: " << __FILE__ << " line:" << __LINE__ << "\n";
    exit(1); 
  }

  typedef Mult<Matrix,PC_Matrix,usePC> MultType;
  typedef typename MultType :: mult_t mult_t;
  // get appropriate mult method 
  mult_t * mult_pc = MultType :: mult_pc;

  typedef double *doubleP;
#ifdef USE_MEMPROVIDER
  int memSize = n*(m+1) 
              + (m*(m+1)/2)
              + n
              + (m+1)
              + m
              + m;

  gmresMem.resize( memSize );
  double *V  = gmresMem.getMem(n*(m+1));
  double *U  = gmresMem.getMem(m*(m+1)/2);
  double *r  = gmresMem.getMem(n);
  double *y  = gmresMem.getMem(m+1);
  double *c  = gmresMem.getMem(m);
  double *s  = gmresMem.getMem(m);
#else 
  double *V  = new double[n*(m+1)];
  double *U  = new double[m*(m+1)/2];
  double *r  = new double[n];
  double *y  = new double[m+1];
  double *c  = new double[m];
  double *s  = new double[m];
#endif
  double **v = new doubleP[m+1];

  // tmp mem for pc mult 
  double * tmp = (usePC) ? (new double[n]) : 0; 

  double error = -1.0;

  for ( int i=0; i<=m; ++i ) 
  {
    v[i]=V+i*n;
  }

  int its=-1;

  {
    double beta, h, rd, dd, nrm2b;
    int j, io, uij, u0j;
    nrm2b = dnrm2(n,b,1);
    
    // global sum 
    nrm2b = comm.sum ( nrm2b );
    
    io=0;
    do  
    { // "aussere Iteration
      ++io;
      mult_pc(A,C,x,r,tmp);
      daxpy(n,-1.,b,1,r,1);
      beta = dnrm2(n,r,1);

      // global sum 
      beta = comm.sum( beta );
      
      dcopy(n,r,1,v[0],1);
      dscal(n,1./beta,v[0],1);

      y[0]=beta;
      j=0;
      uij=0;
      do 
      { // innere Iteration j=0,...,m-1
        u0j=uij;
        mult_pc(A,C,v[j],v[j+1],tmp);
        dgemv(DuneCBlas::Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1);

        // global sum 
        comm.sum( U+u0j, j+1 );
        
        dgemv(DuneCBlas::NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1);

        h = dnrm2(n,v[j+1],1);
        // global sum 
        h = comm.sum( h );
        
        dscal(n,1./h,v[j+1],1);
        
        for ( register int i=0; i<j; ++i ) 
        { // rotiere neue Spalte
          double tmp = c[i]*U[uij]-s[i]*U[uij+1];
          U[uij+1]   = s[i]*U[uij]+c[i]*U[uij+1];
          U[uij]     = tmp;
          ++uij;
        }
        
        { // berechne neue Rotation
          rd     = U[uij];
          dd     = sqrt(rd*rd+h*h);
          c[j]   = rd/dd;
          s[j]   = -h/dd;
          U[uij] = dd;
          ++uij;
        }
        
        { // rotiere rechte Seite y (vorher: y[j+1]=0)
          y[j+1] = s[j]*y[j];
          y[j]   = c[j]*y[j];
        }
        ++j;
        if ( detailed && (comm.rank() == 0))
        {
          std::cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<std::abs(y[j])<<std::endl;
        }
      } 
      while ( j<m && fabs(y[j])>=eps*nrm2b );
      { // minimiere bzgl Y
        dtpsv(UpperTriangle,DuneCBlas::NoTranspose,NotUnitTriangular,j,U,y,1);

        // korrigiere X
        dgemv(DuneCBlas::NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
      }
    } while ( fabs(y[j])>=eps*nrm2b );

    error = std::abs(y[j]);
    
    // R"uckgabe: Zahl der inneren Iterationen
    its = m*(io-1)+j;
  }

  // if right preconditioning then do back solve 
  MultType :: back_solve(n,C,x,tmp);
  
#ifdef USE_MEMPROVIDER
  gmresMem.reset();
#else 
  delete[] V;
  delete[] U;
  delete[] r;
  delete[] y;
  delete[] c;
  delete[] s;
  delete[] v;
#endif

  if( tmp ) 
  {
    delete[] tmp;
  }

  return std::pair<int,double> (its,error);
}

// ============================================================================

template<class CommunicatorType, 
         class Matrix > 
inline 
std::pair<int,double> 
gmres( const CommunicatorType & comm, 
      int m, int n, const Matrix &A, const double *b, double *x, double eps , bool verbose )
{
  return gmres_algo<false> (comm,m,n,A,A,b,x,eps,verbose);
}

template<class CommunicatorType, 
         class Matrix,
         class PC_Matrix > 
inline 
std::pair<int,double> 
gmres( const CommunicatorType & comm, 
      int m, int n, const Matrix &A, const PC_Matrix & C ,
      const double *b, double *x, double eps , bool verbose )
{
  return gmres_algo<true> (comm,m,n,A,C,b,x,eps,verbose);
}

// ============================================================================


#endif // GMRES_BLAS_H
