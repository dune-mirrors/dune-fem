/**************************************************************************
**       Title: ioutils.hh
**    $RCSfile: ioutils.hh,v $
**   $Revision: 1.4 $$Name:  $
**        Date: 28.11.2005
**   Copyright: GPL Bernard Haasdonk
** Description: collection of auxiliary functions for writing 
**              vector/matrix data for further external processing
**************************************************************************/

#ifndef __IOUTILS_HH__
#define __IOUTILS_HH__

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <config.h>

using namespace std;
using namespace Dune;

/*======================================================================*/
/*! 
 *   saveSparseMatrixBinary
 *
 *   simple save routine, which writes a sparse matrix to a binary file. This 
 *   can then be simply read in, e.g. in MATLAB by appropriate reading 
 *   methods. The current version is suited to feop/matrix/spmatrix class
 *   If the old operator/feop/spmatrix class is used, activate code below 
 *   in this file
 * 
 *   param filename the filename to write
 *
 *   param matrix the matrix to write
 *
 *   \return success 0 ( error return 1 is to be implemented)
 */
/*======================================================================*/

// select between implementation for new spmatrix class in 
// operator/matrix/spmatrix or old spmatrix  class in 
// operator/feop/spmatrix

#ifdef USE_OLD_SPARSEMATRIX
// old spmatrix-class:
template <class SparseRowMatrix>
int saveSparseMatrixBinary(const char* filename, SparseRowMatrix& matrix)
{
  // open file for writing  
  ofstream fid(filename, ios::binary | ios::out);    

  // for debugging purposes: write an int and a double
  int magicint =    111;
  double magicdouble = 111.0;
  fid.write((char*)&magicint,sizeof(int));
  fid.write((char*)&magicdouble,sizeof(double));
 
  int nrows = matrix.rows();
  int ncols = matrix.cols();
  int nonzero = matrix.numNonZeros();
      
  // write number of rows and cols and maxnonzeros per row
  fid.write((char*)&nrows,sizeof(int));
  fid.write((char*)&ncols,sizeof(int));
  fid.write((char*)&nonzero,sizeof(int));

  // count nonzero entries
  int totalnonzeros = 0;  
  for (int r=0; r<nrows; r++)
  {
    int nonzero = matrix.numNonZeros(r);    
    for (int c = 0; c < nonzero ; c ++)
    {  
      typename SparseRowMatrix::ColumnIterator it = 
          matrix.rbegin(r);
      for (;it!=matrix.rend(r);++it)
          if (matrix(r,it.col())!=0.0)
              totalnonzeros ++; 
    }
  }
  fid.write((char*)&totalnonzeros,sizeof(int));

  // write all nonzero entries
  for (int r=0; r<nrows; r++)
  {
    int nonzero = matrix.numNonZeros(r);
    for (int c = 0; c < nonzero ; c ++)
    {  
      typename SparseRowMatrix::ColumnIterator it = 
          matrix.rbegin(r);
      for (;it!=matrix.rend(r);++it)
      {
        double entry = matrix(r,it.col());
        int colnum = it.col();
        //! write rownum, colnum and entry
        fid.write((char*)&r,sizeof(int));
        fid.write((char*)&colnum,sizeof(int));
        fid.write((char*)&entry,sizeof(double));        
      }    
    }
  }
  fid.close();  
  return 0; 
}

#else

// new spmatrix-class:
template <class SparseRowMatrix>
int saveSparseMatrixBinary(const char* filename, SparseRowMatrix& matrix)
{
  // open file for writing  
  ofstream fid(filename, ios::binary | ios::out);    
  
  // for debugging purposes and platform check: write an int and a double
  int magicint =    111;
  double magicdouble = 111.0;
  fid.write((char*)&magicint,sizeof(int));
  fid.write((char*)&magicdouble,sizeof(double));
  
  int nrows = matrix.rows();
  int ncols = matrix.cols();
  int nonzero = matrix.numNonZeros();
      
  // write number of rows and cols and maxnonzeros per row
  fid.write((char*)&nrows,sizeof(int));
  fid.write((char*)&ncols,sizeof(int));
  fid.write((char*)&nonzero,sizeof(int));

  // count nonzero entries
  int totalnonzeros = 0;  
  for (int r=0; r<nrows; r++)
  {
    int nonzero = matrix.numNonZeros(r);
    for (int fakeCol = 0; fakeCol < nonzero ; fakeCol ++)
    {  
      int realCol = matrix.realCol(r,fakeCol);
      if ((realCol != SparseRowMatrix::defaultCol) && 
           (matrix(r,realCol)!=0.0))
          totalnonzeros ++; 
    }
  }
  fid.write((char*)&totalnonzeros,sizeof(int));

  // write all nonzero entries
  for (int r=0; r<nrows; r++)
  {
    int nonzero = matrix.numNonZeros(r);
//    for (int c = 0; c < nonzero ; c ++)
//    {  
    
    for (int fakeCol = 0; fakeCol < nonzero ; fakeCol ++)
    {  
      int realCol = matrix.realCol(r,fakeCol);
      if ((realCol != SparseRowMatrix::defaultCol) && 
            (matrix(r,realCol)!=0.0))
      {
        double entry = matrix(r,realCol);
        //! write rownum, colnum and entry
        fid.write((char*)&r,sizeof(int));
        fid.write((char*)&realCol,sizeof(int));
        fid.write((char*)&entry,sizeof(double));        
      }    
    }
  }
  fid.close();  
  return 0; 
}

#endif

/*======================================================================*/
/*! 
 *   saveDofVectorBinary: save dof vector in binary file
 *
 *   \param filename file name to be generated
 *
 *   \param func function to be written
 *
 *   \return 0 success ( error return 1 is to be implemented!)
 */
/*======================================================================*/

template <class DiscreteFunctionType>
int saveDofVectorBinary(const char* filename, DiscreteFunctionType& func)
{  
  typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;    
  
  //const DiscFuncSpaceType& dfsp = func.getFunctionSpace();  

  // open file for writing  
  ofstream fid(filename, ios::binary | ios::out);    

  // for debugging purposes: write an int and a double
  int magicint =    111;
  double magicdouble = 111.0;
  fid.write((char*)&magicint,sizeof(int));
  fid.write((char*)&magicdouble,sizeof(double));
      
  // write number of entries
  int ndofs = func.size();
  fid.write((char*)&ndofs,sizeof(int));
  
  DofIteratorType it = func.dbegin();  
  DofIteratorType eit = func.dend();  
    
  // iterate over all elements defining the function
  for (it = func.dbegin(); it!=eit; ++it)
  {
    double entry = *it;
    fid.write((char*)&entry,sizeof(double));
  } // end element iteration 
  fid.close();
}

#endif
