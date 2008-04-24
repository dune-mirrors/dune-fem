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



namespace Dune
{
/*======================================================================*/
/*! @ingroup DiscFuncIO
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
 *   \param filename the filename to write
 *
 *   \param matrix the matrix to write
 *
 *   \return 1 success or 0 fail
 *
 *   The code for the matlab-reading is realized in load_sparse_matrix_binary.m
 *   within this directory.
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
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    

            // DSM (Dune Sparse Matrix)            
            fid.write("DSM",3);
            
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

            // write end-of-file marker            
            fid.write("EOF",3);

            int status = fid.good();
            
            fid.close();  

            return status; 

          }
    
#else
    
// new spmatrix-class:
    template <class SparseRowMatrix>
    int saveSparseMatrixBinary(const char* filename, SparseRowMatrix& matrix)
          {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    
            
            // write magic number: type of binary file: 
            // DSM (Dune Sparse Matrix)            
            fid.write("DSM",3);

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
            // write end-of-file marker            
            fid.write("EOF",3);

            int status = fid.good();
            
            fid.close();  

            return status; 

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
 *   \return 1 success or 0 fail
 *
 *   The code for the matlab-reading is realized in load_dof_vector_binary.m
 *   within this directory.
 */
/*======================================================================*/
    
    template <class DiscreteFunctionType>
    int saveDofVectorBinary(const char* filename, DiscreteFunctionType& func)
          {  
            typedef typename DiscreteFunctionType::DofIteratorType DofIteratorType;    
            
            //const DiscFuncSpaceType& dfsp = func.getFunctionSpace();  
            
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    

            // write magic number: type of binary file: 
            // DDV (Dune Dof Vector)            
            fid.write("DDV",3);
            
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

            // write end-of-file marker            
            fid.write("EOF",3);

            int status = fid.good();
            
            fid.close();  

            return status; 
            
          }

/*======================================================================*/
/*! 
 *   saveDenseMatrixBinary: save dense matrix in binary file
 *
 *   The matrix class is assumed to give access to the values by 
 *   matrix[nr][nc], values are written as doubles
 *
 *   \param filename file name to be generated
 *
 *   \param matrix matrix to be written
 *
 *   \param nrows number of rows to be written
 *
 *   \param ncols number of columns to be written
 *
 *   \return 1 success or 0 fail
 *
 *   The code for the matlab-reading is realized in load_dune_binary.m
 *   within this directory.
 */
/*======================================================================*/

    template <class DenseRowMatrix>
    int saveDenseMatrixBinary(const char* filename, 
                              DenseRowMatrix& matrix, int nrows, int ncols)
          {
            // open file for writing  
            std::ofstream fid(filename, std::ios::binary | std::ios::out);    
            
            // write magic number: type of binary file: 
            // DDM (Dune Dense Matrix)            
            fid.write("DDM",3);
           
            // for debugging purposes: write an int and a double

            int magicint =    111;
            double magicdouble = 111.0;
            fid.write((char*)&magicint,sizeof(int));
            fid.write((char*)&magicdouble,sizeof(double));
                        
            // write number of rows and cols and maxnonzeros per row
            fid.write((char*)&nrows,sizeof(int));
            fid.write((char*)&ncols,sizeof(int));
            
            // write all entries
            for (int r=0; r<nrows; r++)
            {
              for (int c = 0; c < ncols ; c ++)
              {  
                double entry = matrix[r][c];             
                fid.write((char*)&entry,sizeof(double));        
              }
            }
            
            // write end-of-file marker            
            fid.write("EOF",3);

            int status = fid.good();
            
            fid.close();  

            return status; 
          }
    
  }; // end class MatlabHelper
  
} // namespace Dune

#endif
