namespace Dune 
{

#define EPS 1.0E-15
  
/*****************************/
/*  Constructor(s)           */
/*****************************/
template <class T>
SparseRowMatrix<T>::SparseRowMatrix()
{
  values_ = 0;
  col_ = 0;
  dim_[0] = 0;
  dim_[1] = 0;
  nz_ = 0;
}

template <class T>
SparseRowMatrix<T>::~SparseRowMatrix()
{
  if(values_) delete values_;
  if(col_) delete col_;
}

/***********************************/
/*  Construct from storage vectors */
/***********************************/

template <class T> 
void SparseRowMatrix<T>::
makeSpMat(int rows, int cols, int nz,const T& val )
{
  dim_[0] = rows;
  dim_[1] = cols;

  memSize_ = rows;
  nz_ = nz;

  assert( dim_[0] > 0 );
  assert( dim_[1] > 0 );
  
  values_ = new T [dim_[0]*nz_];
  col_ = new int [dim_[0]*nz_];

  for(int i=0; i<dim_[0]*nz_; i++)
  { 
    values_[i] = val;
    col_[i] = -1;
  }
}

template <class T> 
SparseRowMatrix<T>::
SparseRowMatrix(int rows, int cols, int nz, const T& val)
{
  makeSpMat(rows,cols,nz,val);
}

template <class T> 
T SparseRowMatrix<T>::operator()(int row, int col) const
{
  for (int i=0; i<nz_; i++)
  { 
    if(col_[row*nz_ +i] == col)
    {
      return values_[row*nz_ +i];
    }
  }
  return 0;
}

template <class T> 
int SparseRowMatrix<T>::colIndex(int row, int col)
{
  int whichCol = -1;
  int thisCol = 0;
  for(int i=0; i<nz_; i++)
  {
    thisCol = col_[row*nz_ +i];
    if(col == thisCol) return i;
    if(thisCol == -1) return i;
  }
  
  //assert(whichCol != -1);
  if(whichCol < 0) 
  {
    for(int i=0; i<nz_; ++i) 
    {
      std::cout << col_[row*nz_ +i] << " ";
    }
    std::cout << std::endl;
  }
  return whichCol;
}

template <class T> 
bool SparseRowMatrix<T>::find (int row, int col) const
{
  int thisCol = 0;
  for(int i=0; i<nz_; ++i)
  {
    thisCol = col_[row*nz_ +i];
    if(col == thisCol) return true;
    if(thisCol == -1) return false;
  }
  return false;
}

template <class T> 
void SparseRowMatrix<T>::clear()
{
  for(register int i=0; i<dim_[0]*nz_; ++i)
  {
    values_ [i] = 0;
    col_[i] = -1;
  }
}

template <class T> 
void SparseRowMatrix<T>::set(int row, int col, T val)
{
  int whichCol = colIndex(row,col);
  assert( whichCol >= 0);

  {
    values_[row*nz_ + whichCol] = val; 
    col_[row*nz_ + whichCol] = col;
  }
}

template <class T> 
void SparseRowMatrix<T>::add(int row, int col, T val)
{
  int whichCol = colIndex(row,col);
  assert(whichCol >= 0);
  values_[row*nz_ + whichCol] += val; 
  col_[row*nz_ + whichCol] = col;
}

template <class T> 
void SparseRowMatrix<T>::multScalar(int row, int col, T val)
{
  int whichCol = colIndex(row,col);
  if(whichCol < 0)
  {
    std::cout << " error \n";
  }
  else
  {
    values_[row*nz_ + whichCol] *= val; 
    col_[row*nz_ + whichCol] = col;
  }
}
/***************************************/
/*  Matrix-MV_Vector multiplication    */
/***************************************/
template <class T> template <class VECtype>
void SparseRowMatrix<T>::mult(const VECtype *x, VECtype *ret) const
{
  multOEM(x,ret);
}

template <class T> template <class VECtype>
void SparseRowMatrix<T>::multOEM(const VECtype *x, VECtype *ret) const
{
  for(register int row=0; row<dim_[0]; ++row)
  {
    T sum = 0;
    int thisCol = row*nz_ ;
    const T * localValues = &values_[thisCol];
    for(int col=0; col<nz_; ++col)
    {
      int realCol = col_[ thisCol ];
      if ( realCol < 0 ) break;
      sum += localValues[col] * x[ realCol ];
      ++thisCol; 
    }
    ret[row] = sum;
  }
  return; 
}

template <class T> template <class VECtype>
void SparseRowMatrix<T>::multOEM_t(const VECtype *x, VECtype *ret) const
{
  for(register int col=0; col<dim_[1]; ++col)
  {
    ret[col] = 0.0;
  }

  for(register int row=0; row<dim_[0]; ++row)
  {
    for(register int col=0; col<nz_; ++col)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[ thisCol ];
      if ( realCol < 0 ) continue;
      ret[realCol] += values_[thisCol] * x[ row ];
    }
  }
  return; 
}

/***************************************/
/*  Matrix-MV_Vector multiplication    */
/***************************************/
#if 0
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::apply(const DiscFuncType &f, DiscFuncType &ret) const 
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType ret_it = ret.dbegin(); 
  const DofIteratorType f_it = f.dbegin(); 

  for(int row=0; row<dim_[0]; row++)
  {
    (*ret_it) = 0.0;
    
    //! DofIteratorType schould be the same 
    for(int col=0; col<nz_; col++)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[thisCol];
      
      if( realCol < 0 ) continue;        
      (*ret_it) += values_[thisCol] * (f_it[realCol]);
    }

    ++ret_it;
  } 

  return; 
}
#endif

template <class T> template <class DiscFType , class DiscFuncType>
void SparseRowMatrix<T>::apply(const DiscFType &f, DiscFuncType &ret) const 
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  
  typedef typename DiscFuncType::ConstDofIteratorType ConstDofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType ret_it = ret.dbegin(); 
  ConstDofIteratorType f_it = f.dbegin(); 

  for(int row=0; row<dim_[0]; row++)
  {
    (*ret_it) = 0.0;
    
    //! DofIteratorType schould be the same 
    for(int col=0; col<nz_; col++)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[thisCol];
      
      if( realCol < 0 ) continue;        
      (*ret_it) += values_[thisCol] * (f_it[realCol]);
    }

    ++ret_it;
  } 

  return; 
}


// apply to tranpose matrix 
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::apply_t(const DiscFuncType &f, DiscFuncType &ret) const 
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  
  int level = f.getFunctionSpace().getGrid().maxlevel();

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType ret_it = ret.dbegin(); 
  const DofIteratorType f_it = f.dbegin(); 

  for(int row=0; row<dim_[0]; row++)
  {
    (*ret_it) = 0.0;
    
    //! DofIteratorType schould be the same 
    for(int col=0; col<nz_; col++)
    {
      int thisCol = col * nz_ + row;
      int realCol = col_[thisCol];
      
      if( realCol < 0 ) continue;        
      (*ret_it) += values_[thisCol] * (f_it[realCol]);
    }

    ++ret_it;
  } 

  return; 
}


template <class T> 
void SparseRowMatrix<T>::print(std::ostream& s) const
{
  for(int row=0; row<dim_[0]; row++)
  {
    for(int col=0; col<dim_[1]; col++)
    {
      s << (*this)(row,col) << " ";
    }
    s << "\n";
  }
  return; 
}

template <class T> 
void SparseRowMatrix<T>::printReal(std::ostream& s) const
{
  for(int row=0; row<dim_[0]; row++)
  {
    for(int col=0; col<nz_; col++)
    {
      s << values_[row*nz_ + col] << " ";
    }
    s << "\n";
  }
  return; 
}

template <class T> 
void SparseRowMatrix<T>::kroneckerKill(int row, int col) 
{
  unitRow(row);
  unitCol(col);
} 

template <class T> 
void SparseRowMatrix<T>::unitRow(int row) 
{
  for(int i=1; i<nz_; i++) 
  {
    values_[row*nz_ + i] = 0.0;
    col_[row*nz_ + i] = -1; 
  }
  values_[row*nz_] = 1.0;
  col_[row*nz_] = row;
} 

template <class T> 
void SparseRowMatrix<T>::unitCol(int col) 
{
  for(int i=0; i<dim_[0]; i++) 
    if (i != col) set(i,col,0.0); 
    else set(col,col,1.0); 
} 

template <class T> 
void SparseRowMatrix<T>::checkSym() 
{
  double val;
  for(int i=0; i<this->size(0); i++)
  {
    for(int j=0; j<this->size(0); j++)
    {
      val = this->operator() (i,j);
      if(std::abs(val - this->operator() (j,i)) > 1E-10)
      {
        std::cout << val << " " << this->operator() (j,i) << " val \n";
      }
    }  
  }
} 


// diagonal conditioning  
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::diagCond(DiscFuncType &rhs)  
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType it = rhs.dbegin(); 

  for(int row=0; row<dim_[0]; row++)
  {
    double diag_1 = this->operator() (row,row);
    if(std::abs(diag_1) > 0.0)
    {
      diag_1 = 1.0/diag_1;
    }
    else 
    {
      diag_1 = 1.0;
    }

    // multipy rhs with 1/diag 
    (*it) *= diag_1;
    ++it;
   
    //! DofIteratorType schould be the same 
    for(int col=0; col<nz_; col++)
    {
      int thisCol = nz_ * row + col;
      //int realCol = col_[thisCol];
      //if( realCol < 0 ) continue;        
      values_[thisCol] *= diag_1;
    }
  } 

  return; 

}

// diagonal conditioning  
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::getDiag(const ThisType & mass,
                                 const ThisType & B, 
                                 DiscFuncType &diag) const  
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType it = diag.dbegin(); 

  assert( this->size(0) == B.size(1) );
  assert( mass.size(0) == mass.size(1) ); 
  assert( mass.size(0) == this->size(1) );

  for(int row=0; row<this->size(0); row++)
  {
    T sum = 0.0;
    for(register int col=0; col<nz_; ++col)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[ thisCol ];
      if ( realCol < 0 ) continue;
      
      double diag = mass(realCol,realCol);
      sum += diag * values_[ thisCol ] * B(realCol,row);
    }

    (*it) = sum;
    ++it;
  }
  return; 
}

// diagonal conditioning  
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::getDiag(const ThisType & B, 
                                 DiscFuncType &diag) const  
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType it = diag.dbegin(); 

  assert( this->size(0) == B.size(1) );

  for(int row=0; row<this->size(0); row++)
  {
    T sum = 0.0;
    for(register int col=0; col<nz_; ++col)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[ thisCol ];
      if ( realCol < 0 ) continue;
      sum += values_[ thisCol ] * B(realCol,row);
    }

    (*it) = sum;
    ++it;
  }
  return; 
}

// diagonal conditioning  
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::getDiag(DiscFuncType &diag) const  
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType it = diag.dbegin(); 

  for(int row=0; row<this->size(0); row++)
  {
    (*it) = (*this)(row,row); 
    ++it;
  }
  return; 
}

// add diagonal to given discrete function   
template <class T> template <class DiscFuncType>
void SparseRowMatrix<T>::addDiag(DiscFuncType &diag) const  
{
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType it = diag.dbegin(); 

  for(int row=0; row<this->size(0); row++)
  {
    (*it) += (*this)(row,row); 
    ++it;
  }
  return; 
}
// resize with rows = cols = newSize  
template <class T>
void SparseRowMatrix<T>::resize (int newSize)  
{
  resize(newSize,newSize);
}

// resize matrix 
template <class T>
void SparseRowMatrix<T>::resize (int newRow, int newCol)  
{
  if(newRow != this->size(0) && newCol != this->size(1))
  {
    if(newRow != memSize_)
    {
      if(values_) delete [] values_;
      if(col_) delete [] col_;

      T tmp = 0;
      makeSpMat(newRow,newCol,nz_,tmp);
    }
    else 
    {
      assert(newRow > 0);
      dim_[0] = newRow;
      dim_[1] = newCol;
    }
  }
}

template <class T> 
void SparseRowMatrix<T>::multiply(const SparseRowMatrix<T> & B,
    SparseRowMatrix<T> & res) const 
{
  //res.resize( B.size(0) );
  res.clear();
  //assert( res.NumNonZeros() == B.NumNonZeros() );

  //std::cout << res.NumNonZeros() << "\n";
  for(int row=0; row<this->size(0); row++)
  {
    for(int col=0; col<B.size(1); col++)
    {
      T sum = 0;
      for(int k=0; k<B.size(0); k++)
      {
        sum += (*this)(row,k) * B(k,col);
      }

      if( std::abs(sum) > 0.0)
      {
        //std::cout << "add ("<<row<<","<<col<<")\n";
        res.add(row,col,sum);
      }
    }
  }

  //res.print(std::cout);
}

template <class T> 
void SparseRowMatrix<T>::add(const SparseRowMatrix<T> & B)
{
  assert( this->size(0) == B.size(0) );
  assert( nz_ == B.nz_ );
  
  for(register int i=0; i<dim_[0]*nz_; ++i)
  {
    assert( col_ [i] == B.col_ [i] );
    values_ [i] += B.values_[i];
  }
}

} // end namespace Dune
