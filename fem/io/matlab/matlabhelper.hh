/**************************************************************************
**       Title: matlabhelper.hh
**    $RCSfile: matlabhelper.hh,v $
**   $Revision: 1.4 $$Name:  $
**        Date: 28.11.2005
**   Copyright: GPL Bernard Haasdonk
** Description: collection of auxiliary functions for writing 
**              vector/matrix data for further external processing, e.g. in
**              MATLAB
**************************************************************************/

#ifndef __MATLABHELPER_HH__
#define __MATLABHELPER_HH__

#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <config.h>

using namespace std;
using namespace Dune;


namespace Dune
{
/*======================================================================*/
/*!
 *  \class MatlabHelper
 *  \brief The MatlabHelper class provides functionality for exporting
 *         Dune Structures to Matlab.
 *
 *   Currently only two methods are provides, which perform a binary 
 *   saving of SparseRowMatrix and DiscreteFunction's DOFVectors to a
 *   binary format, which can be read by appropriate MATLAB-files
 */
/*======================================================================*/


  class MatlabHelper
  {
  public:
    
/*======================================================================*/
/*! 
 *   saveSparseMatrixBinary
 *
 *   simple save routine, which writes a sparse matrix to a binary file. This 
 *   can then be simply read in, e.g. in MATLAB by appropriate reading 
 *   methods. The current version is suited to feop/matrix/spmatrix class
 *   The MATLAB-code for the reading method is the function
 *   load_sparse_matrix.m 
 * 
 *   param filename the filename to write
 *
 *   param matrix the matrix to write
 *
 *   \return success 0 ( error return 1 is to be implemented)
 */
/*======================================================================*/

/*======================================================================*/
/* MATLAB-Code for reading: save the following as load_sparse_matrix.m
function A = load_sparse_matrix(filename)
%function A = load_sparse_matrix(filename)
%
% load binary file, which represents a sparse matrix, e.g. generated
% by saveSparseMatrix
%
% format: 
% magic numbers: int 111, double 111
% number of rows, cols and maxnoonzero_per_row
% number of total_nonzeros
% for 0.. total_nonzeros-1 : triples (int r,int c,double v)   
% where r,c start from 0

% Bernard Haasdonk 15.12.2006

  fid = fopen(filename,'r');
  %fid = fopen(filename,'r','ieee-be');
  
  magicint = fread(fid,1,'int');
  magicdouble = fread(fid,1,'double');
  
  if (magicint~=111) | (magicdouble~=111.0)
    error('magic numbers not read correctly!');
  end;
  
  nrows = fread(fid,1,'int');
  ncols = fread(fid,1,'int');
  nnonzeros = fread(fid,1,'int');
  ntotalnonzeros = fread(fid,1,'int');
  
  disp(['generating ',num2str(nrows),'x',num2str(ncols),...
	' sparse matrix with ',num2str(ntotalnonzeros),' totalnonzeros.']);
  A = sparse(nrows,ncols,nnonzeros);
  
  for i=1:ntotalnonzeros
    row = fread(fid,1,'int');
    col = fread(fid,1,'int');
    val = fread(fid,1,'double');
    A(row+1,col+1) = val;
  end;  
*/ /* END OF MATLAB CODE*/
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
 *
 * The code for the Matlab-reading method can be extracted from this file.
 */
/*======================================================================*/

/*======================================================================*/
/* MATLAB-Code for reading: save the following as load_dof_vector.m

function v = load_dof_vector(filename)
%function A = load_dof_vector(filename)
%
% load binary file, which represents a vector of double values
%
% format: 
% magic numbers: int 111, double 111
% number of entries
% for 0.. numer_of_entries-1 double v   

% Bernard Haasdonk 15.12.2006

  fid = fopen(filename,'r');
  %fid = fopen(filename,'r','ieee-be');
  
  magicint = fread(fid,1,'int');
  magicdouble = fread(fid,1,'double');
  
  if (magicint~=111) | (magicdouble~=111.0)
    error('magic numbers not read correctly!');
  end;
  
  nentries = fread(fid,1,'int');
  
  disp(['generating and reading vector with ',num2str(nentries),' entries.']);

  v = fread(fid,nentries,'double');
*/ /* END OF MATLAB CODE*/
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
    
  }; // end class MatlabHelper
  
} // namespace Dune

#endif
