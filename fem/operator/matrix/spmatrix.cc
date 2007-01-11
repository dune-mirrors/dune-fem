namespace Dune 
{

#define EPS 1.0E-15
  
/*****************************/
/*  Constructor(s)           */
/*****************************/
template <class T>
SparseRowMatrix<T>::SparseRowMatrix(double omega) : omega_(omega)
{
  values_ = 0;
  col_ = 0;
  dim_[0] = 0;
  dim_[1] = 0;
  memSize_ = 0;
  nz_ = 0;
  nonZeros_ = 0;
}

template <class T>
void SparseRowMatrix<T>::removeObj()
{
  if(values_) delete [] values_;
  if(col_) delete [] col_;
  if(nonZeros_) delete [] nonZeros_;
}

template <class T>
SparseRowMatrix<T>::~SparseRowMatrix()
{
  removeObj();
}

/***********************************/
/*  Construct from storage vectors */
/***********************************/
template <class T> 
void SparseRowMatrix<T>::
reserve(int rows, int cols, int nz,const T& val )
{
  if( (rows == dim_[0]) && (cols == dim_[1]) && (nz == nz_)) 
  {
    clear();
    return;
  }

  removeObj();

  values_ = new T [ rows*nz ];
  col_    = new int [ rows*nz ];
  nonZeros_ = new int [ rows ];

  assert( values_ );
  assert( col_ );
  assert( nonZeros_ );

  dim_[0] = rows;
  dim_[1] = cols;

  memSize_ = rows;
  nz_ = nz;
  // add first col for offset
  nz_ += firstCol ;

  assert( dim_[0] > 0 );
  assert( dim_[1] > 0 );
  
  // make resize 
  newValues_.resize( nz_ );
  
  // only reserve for indices 
  newIndices_.reserve( nz_ );

  // set all values to default values 
  clear();
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
  if(newRow != this->size(0))
  {
    int memHalf = (int) memSize_/2;
    if((newRow > memSize_) || (newRow < memHalf))
    {
      T tmp = 0;
      T * oldValues = values_;       values_ = 0;
      int * oldCol  = col_;          col_ = 0;
      int * oldNonZeros = nonZeros_; nonZeros_ = 0;
      int oldNz = nz_;
      int copySize = std::min( dim_[0] , newRow );
      int oldSize = dim_[0];

      // reserve new memory 
      reserve(newRow,newCol,nz_,tmp);

      if( (oldSize > 0) && (oldNz > 0 ))
      {
        std::memcpy(values_  , oldValues  , copySize * nz_ * sizeof(T) );
        std::memcpy(col_     , oldCol     , copySize * nz_ * sizeof(int) );
        std::memcpy(nonZeros_, oldNonZeros, copySize * sizeof(int) );
      }

      delete [] oldValues;
      delete [] oldCol; 
      delete [] oldNonZeros; 
    }
    else 
    {
      assert(newRow > 0);
      dim_[0] = newRow;
      dim_[1] = newCol;
    }
  }

  assert( this->size(0)  == newRow );
  assert( this->size(1)  == newCol );
}

template <class T> 
SparseRowMatrix<T>::
SparseRowMatrix(int rows, int cols, int nz, const T& val)
{
  reserve(rows,cols,nz,val);
}

template <class T> 
T SparseRowMatrix<T>::operator()(int row, int col) const
{
  assert( row >= 0 );
  assert( (row < dim_[0]) ? 1 : (std::cout << row << " bigger " << dim_[0] <<"\n", 0));

  const int nonZ = nonZeros_[row]; 
  int thisCol = row*nz_;
  for (int i=firstCol; i<nonZ; ++i)
  { 
    if(col_[thisCol] == col)
    {
      return values_[thisCol];
    }
    ++thisCol; 
  }
  return 0;
}

template <class T> 
int SparseRowMatrix<T>::colIndex(int row, int col)
{
  assert( row >= 0 );
  assert( row < dim_[0] );

  int whichCol = defaultCol;
  int thisCol = 0;
  for(int i=firstCol; i<nz_; ++i)
  {
    thisCol = col_[row*nz_ +i];
    if(col == thisCol) return i;
    if(thisCol == defaultCol ) 
    {
      ++nonZeros_[row];
      return i;
    }
  }
  
  if(whichCol == defaultCol ) 
  {
    std::cout << "Writing colIndex for, nz = " << nz_ <<  " , " << col << "\n";
    for(int i=firstCol; i<nz_; ++i) 
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
  for(int i=firstCol; i<nz_; ++i)
  {
    thisCol = col_[row*nz_ +i];
    if(col == thisCol) return true;
    if(thisCol == defaultCol ) return false;
  }
  return false;
}

template <class T> 
void SparseRowMatrix<T>::clear()
{
  T init = 0;
  for(register int i=0; i<dim_[0]*nz_; ++i)
  {
    values_ [i] = init;
    col_[i] = defaultCol;
  }

  for(register int i=0; i<dim_[0]; ++i)
  {
    nonZeros_[i] = 0;
  }
}

template <class T> 
void SparseRowMatrix<T>::clearRow(int row)
{
  assert( nonZeros_ );
  assert( values_ );
  assert( col_ );

  nonZeros_[row] = firstCol;
  
  int col = row * nz_;
  for(int i=0; i<nz_; ++i)
  {
    values_ [col] = 0;
    col_[col] = defaultCol;
    ++col;
  }
}

template <class T> 
void SparseRowMatrix<T>::resort()
{
  const int nRows = rows();
  for(int row=0; row<nRows; ++row)
  {
    resortRow(row);
  }
}

template <class T> 
void SparseRowMatrix<T>::resortRow(const int row)
{
  newIndices_.resize(0);
  int thisCol = row * nz_;
  
  for(int col=0; col<nz_; ++col)
  {
    int realCol =  col_[ thisCol + col ] ;
    if( realCol > defaultCol )
    {
      newIndices_.push_back( realCol );
    }
  }
   
  // set number of non zeros for row 
  const int nZero = newIndices_.size();

  // nonZeros should be already at right size 
  assert( nonZeros_[row] == nZero );
  //nonZeros_[row] = nZero;
  //std::cout << "found nz = " << nZero << "\n";

  // make values cache efficient 
  std::sort( newIndices_.begin(), newIndices_.end() );
  for(int col=0; col<nZero; ++col)
  {
    int val = col_[ thisCol + col ];
    T value = values_[ thisCol + col ];
    for(int j=0; j<nZero; ++j) 
    {
      if( newIndices_[j] == val ) 
        newValues_[j] = value; 
    }
  }

  for(int col=0; col<nZero; ++col)
  {
    values_[ thisCol ] = newValues_[col];
    col_[ thisCol ] = newIndices_[col];
    ++thisCol;
  }
}

template <class T> 
void SparseRowMatrix<T>::set(int row, int col, T val)
{
  int whichCol = colIndex(row,col);
  assert( whichCol != defaultCol );

  {
    values_[row*nz_ + whichCol] = val; 
    col_[row*nz_ + whichCol] = col;
  }
}

template <class T> 
void SparseRowMatrix<T>::add(int row, int col, T val)
{
  int whichCol = colIndex(row,col);
  assert( whichCol != defaultCol );
  values_[row*nz_ + whichCol] += val; 
  col_[row*nz_ + whichCol] = col;
}

template <class T> 
void SparseRowMatrix<T>::multScalar(int row, int col, T val)
{
  int whichCol = colIndex(row,col);
  
  if(whichCol == defaultCol )
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
T SparseRowMatrix<T>::multOEMRow(const VECtype *x, const int row) const
{
  T sum = 0;
  int thisCol = row*nz_ + firstCol ;
  const T * localValues = &values_[thisCol];
  const int nonZero = nonZeros_[row];
  for(int col = firstCol ; col<nonZero; ++col)
  {
    int realCol = col_[ thisCol ];
    assert( realCol > defaultCol );
    sum += localValues[col] * x[ realCol ];
    ++thisCol; 
  }
  return sum; 
}

template <class T> template <class VECtype>
void SparseRowMatrix<T>::multOEM(const VECtype *x, VECtype *ret) const
{
  for(register int row=0; row<dim_[0]; ++row)
  {
    ret[row] = multOEMRow( x, row ); 
  }
  return; 
}

template <class T> template <class VECtype>
void SparseRowMatrix<T>::multOEMAdd(const VECtype *x, VECtype *ret) const
{
  for(register int row=0; row<dim_[0]; ++row)
  {
    ret[row] += multOEMRow( x, row ); 
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
    const int nonZero = nonZeros_[row];
    for(register int col=0; col<nonZero; ++col)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[ thisCol ];
      assert( realCol > defaultCol );
      ret[realCol] += values_[thisCol] * x[ row ];
    }
  }
  return; 
}

/***************************************/
/*  Matrix-MV_Vector multiplication    */
/***************************************/

template <class T> template <class DiscFType , class DiscFuncType>
void SparseRowMatrix<T>::apply(const DiscFType &f, DiscFuncType &ret) const 
{
  multOEM(f.leakPointer(),ret.leakPointer());
  /*
  
  typedef typename DiscFuncType::DofIteratorType DofIteratorType;  
  typedef typename DiscFuncType::ConstDofIteratorType ConstDofIteratorType;  

  //! we assume that the dimension of the functionspace of f is the same as
  //! the size of the matrix 
  DofIteratorType ret_it = ret.dbegin(); 
  ConstDofIteratorType f_it = f.dbegin(); 

  for(int row=0; row<dim_[0]; row++)
  {
    (*ret_it) = 0.0;
    
    T sum = 0;
    int thisCol = row*nz_ + firstCol ;

    //! DofIteratorType schould be the same 
    const T * localValues = &values_[thisCol];
    const int nonZero = nonZeros_[row];
    for(int col=firstCol; col<nz_; col++)
    {
      int thisCol = row*nz_ + col;
      int realCol = col_[thisCol];
      
      if( realCol == defaultCol ) continue;        
      (*ret_it) += values_[thisCol] * (f_it[realCol]);
    }

    ++ret_it;
  } 
  */

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
    for(int col=firstCol; col<nz_; col++)
    {
      int thisCol = col * nz_ + row;
      int realCol = col_[thisCol];
      
      if( realCol == defaultCol ) continue;        
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
    col_[row*nz_ + i] = defaultCol; 
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

  for(int row=0; row<this->size(0); ++row)
  {
    T sum = 0.0;
    int thisCol = row*nz_;
    const T * localValues = &values_[thisCol];
    const int nZeros = nonZeros_[row];
    for(register int col=0; col<nZeros; ++col)
    {
      int realCol = col_[ thisCol ];
      assert( realCol != defaultCol );
      double diag = mass(realCol,realCol);
      sum += diag * localValues[col] * B(realCol,row);
      ++thisCol;
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

template <class T> 
void SparseRowMatrix<T>::ssorPrecondition(const T* u, T* x) const 
{
  const double omega = omega_;

  // (D - omega E) x = x_old (=u)  
  for(int row=0; row<dim_[0]; ++row)
  {
    double diag=1.0, dot=0.0;
    // get row stuff 
    int thisCol = row*nz_ + firstCol ;
    const T * localValues = &values_[thisCol];
    const int nonZero = nonZeros_[row];
    for(int col = firstCol ; col<nonZero; ++col)
    {
      const int realCol = col_[ thisCol ];
      assert( realCol > defaultCol );
      
      if (realCol < row) 
      {
        dot += localValues[col] * x[realCol];
      }
      else if (realCol == row) 
      {
        diag = localValues[col];
      }
      ++thisCol;
    }

    x[row] = (u[row] - omega*dot) / diag;
  }

  // D^{-1} (D - omega F) x = x_old (=x)
  for(int row=dim_[0]-1; row>=0; --row)
  {
    double diag=1.0, dot=0.0;
    int thisCol = row*nz_ + firstCol ;
    const T * localValues = &values_[thisCol];
    const int nonZero = nonZeros_[row];
    for(int col = firstCol ; col<nonZero; ++col)
    {
      const int realCol = col_[ thisCol ];
      assert( realCol > defaultCol );
      
      if (realCol > row) 
      {
        dot += localValues[col] * x[realCol];
      }
      else if (realCol == row) 
      {
        diag = localValues[col];
      }
      ++thisCol;
    }
    x[row] -= omega * dot / diag;
  }
}



} // end namespace Dune
