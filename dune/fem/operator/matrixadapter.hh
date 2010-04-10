/**************************************************************************
**       Title: Collection of simple matrix-adapter classes
**    $RCSfile$
**   $Revision$$Name$
**       $Date$
**   Copyright: GPL $Author$
** Description: matrix adapters for 
**              FieldMatrix ==> LocalMatrix-functionality used in FEOp
**              .....
**
**************************************************************************/

#ifndef DUNE_MATRIXADAPTER_HH
#define DUNE_MATRIXADAPTER_HH

#include <dune/fem/version.hh>

namespace Dune
{
  
/*======================================================================*/
/*!
 *  \class FieldMatrixAdapter
 *  \brief Extend Fieldmatrix class to have certain functions as required 
 *         by FEOp, etc. 
 *
 *  A FieldMatrixAdapter provides a  
 *     constructor without arguments
 *     rows(), cols() methods
 *     add(rown, coln, value) allows writable additive access to 
 *         ij-th component.
 *     clear() method
 *     const double& operator() (int i, int j) const : element value 
 *                  readable access
 */
/*======================================================================*/

template <class FieldMatrixImp>
class FieldMatrixAdapter
{
public:  
  typedef FieldMatrixImp FieldMatrixType;
  typedef typename FieldMatrixType::Iterator Iterator;
  
  //! constructor
  FieldMatrixAdapter() DUNE_VERSION_DEPRECATED(1,2,remove)
          : mat_()
        {};

  //! destructor    
  ~FieldMatrixAdapter()
        {};  

  //! determine number of rows
  inline int rows()
        {
          return mat_.rows;
        }
  
  //! determine number of columns
  inline int cols()
        {
          return mat_.cols;
        }
  
  //! addition of element
  inline void add(int nrow, int ncol, double value)
        {
          mat_[nrow][ncol]+= value;
        }  

  //! clear matrix: set to zero
  inline void clear()
        {
          Iterator it = mat_.begin();
          for (;it!=mat_.end();++it)
              *it=0.0;
        }
  
  //!  Const index operator. 
  inline const double& operator() (int i, int j) const
        {
          return mat_[i][j];
        }
  
  
private:
  FieldMatrixType mat_;
};
 
}; // end namespace Dune

#endif



