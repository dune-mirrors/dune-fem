//  -*- C++ -*-

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

template< class Matrix >
inline 
std::pair<int,double> 
gmres( int m, int N, const Matrix &A, const double *b, double *x, double eps );


template< class Matrix >
inline 
std::pair<int,double> 
gmres( int m, int N, const Matrix &A, const double *b, double *x, double eps,
       bool detailed );


// ============================================================================


#include "cblas.h"


template< class Matrix >
inline 
std::pair<int,double> 
gmres( int m, int n, const Matrix &A, const double *b, double *x, double eps,
       bool detailed ) 
{
  if ( n<=0 ) 
  {
    std::cerr << "WARNING: n = " << n << " in gmres, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  typedef double *doubleP;
  double *V  = new double[n*(m+1)];
  double *U  = new double[m*(m+1)/2];
  double *r  = new double[n];
  double *y  = new double[m+1];
  double *c  = new double[m];
  double *s  = new double[m];
  double **v = new doubleP[m+1];

  double error = -1.0;

  for ( int i=0; i<=m; ++i ) v[i]=V+i*n;
  int its=-1;
  {
    double beta, h, rd, dd, nrm2b;
    int j, io, uij, u0j;
    nrm2b=dnrm2(n,b,1);
    
    io=0;
    do  
    { // "aussere Iteration
      ++io;
      mult(A,x,r);
      daxpy(n,-1.,b,1,r,1);
      beta=dnrm2(n,r,1);
      dcopy(n,r,1,v[0],1);
      dscal(n,1./beta,v[0],1);

      y[0]=beta;
      j=0;
      uij=0;
      do 
      { // innere Iteration j=0,...,m-1
        u0j=uij;
        mult(A,v[j],v[j+1]);
        dgemv(Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1);
        dgemv(NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1);
        h=dnrm2(n,v[j+1],1);
        dscal(n,1./h,v[j+1],1);
        for ( int i=0; i<j; ++i ) 
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
        
        if ( detailed )
        {
          std::cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<std::abs(y[j])<<std::endl;
        }
        
      } while ( j<m && fabs(y[j])>=eps*nrm2b );
      
      { // minimiere bzgl Y
        dtpsv(UpperTriangle,NoTranspose,NotUnitTriangular,j,U,y,1);
        // korrigiere X
        dgemv(NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
      }
    } while ( fabs(y[j])>=eps*nrm2b );

    error = std::abs(y[j]);
    // R"uckgabe: Zahl der inneren Iterationen
    its = m*(io-1)+j;
  }

  delete[] V;
  delete[] U;
  delete[] r;
  delete[] y;
  delete[] c;
  delete[] s;
  delete[] v;

  return std::pair<int,double> (its,error);
}

template< class Matrix , class PC_Matrix >
inline
std::pair<int,double> 
gmres_pc ( int m, int n, const Matrix &A, const PC_Matrix & C, const double *rhs , double *x, double eps,
       bool detailed ) 
{
  if ( n<=0 )
  {
    std::cerr << "WARNING: n = " << n << " in gmres_pc, file: " << __FILE__ << " line:" << __LINE__ << "\n";
    return std::pair<int,double> (-1,0.0);
  }

  double * newRhs = new double[n];
  double * tmp = new double[n];
  for(register int i=0; i<n; ++i) newRhs[i] = 0.0;

  mult(C,rhs,newRhs);
  const double * b = newRhs;
  
  typedef double *doubleP;
  double *V  = new double[n*(m+1)];
  double *U  = new double[m*(m+1)/2];
  double *r  = new double[n];
  double *y  = new double[m+1];
  double *c  = new double[m];
  double *s  = new double[m];
  double **v = new doubleP[m+1];

  double error = -1.0;

  for ( int i=0; i<=m; ++i ) 
  {
    v[i]=V+i*n;
  }

  int its=-1;

  {
    double beta, h, rd, dd, nrm2b;
    int j, io, uij, u0j;
    nrm2b=dnrm2(n,b,1);
    
    io=0;
    do  
    { // "aussere Iteration
      ++io;
      for(register int k=0; k<n; ++k) tmp[k] = 0.0;
      mult_pc(A,C,x,r,tmp);
      daxpy(n,-1.,b,1,r,1);
      beta=dnrm2(n,r,1);
      dcopy(n,r,1,v[0],1);
      dscal(n,1./beta,v[0],1);

      y[0]=beta;
      j=0;
      uij=0;
      do 
      { // innere Iteration j=0,...,m-1
        u0j=uij;
        for(register int k=0; k<n; ++k) tmp[k] = 0.0;
        mult_pc(A,C,v[j],v[j+1],tmp);
        dgemv(Transpose,n,j+1,1.,V,n,v[j+1],1,0.,U+u0j,1);
        dgemv(NoTranspose,n,j+1,-1.,V,n,U+u0j,1,1.,v[j+1],1);
        h=dnrm2(n,v[j+1],1);
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
        if ( detailed )
        {
          std::cout<<"gmres("<<m<<")\t"<<io<<"\t"<<j<<"\t"<<std::abs(y[j])<<std::endl;
        }
      } 
      while ( j<m && fabs(y[j])>=eps*nrm2b );
      { // minimiere bzgl Y
        dtpsv(UpperTriangle,NoTranspose,NotUnitTriangular,j,U,y,1);
        // korrigiere X
        dgemv(NoTranspose,n,j,-1.,V,n,y,1,1.,x,1);
      }
    } while ( fabs(y[j])>=eps*nrm2b );

    error = std::abs(y[j]);
    
    // R"uckgabe: Zahl der inneren Iterationen
    its = m*(io-1)+j;
  }

  delete[] V;
  delete[] U;
  delete[] r;
  delete[] y;
  delete[] c;
  delete[] s;
  delete[] v;
  delete[] newRhs;
  delete[] tmp;
  
  return std::pair<int,double> (its,error);
}


// ============================================================================

template< class Matrix >
inline 
std::pair<int,double> 
gmres( int m, int n, const Matrix &A, const double *b, double *x, double eps ){
  return gmres(m,n,A,b,x,eps,false);
}

// ============================================================================


#endif // GMRES_BLAS_H
