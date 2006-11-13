#ifndef DUNE_BLOCKSPMATRIX_HH
#define DUNE_BLOCKSPMATRIX_HH

//- system includes 
#include <stack>

//- local includes 
#include "spmatrix.hh"

namespace Dune
{

//*****************************************************************
//
//  --BlockMatrix
//  
//! Compressed row sparse matrix, where only the nonzeros of a row a 
//! keeped
//*****************************************************************
template <class T>
class DenseMatrix 
{
public: 
  typedef T Ttype;  //! remember the value type

private:
  typedef std::vector < std::vector < T > > MatrixType; 
  typedef  std::vector < T > RowType;

  MatrixType matrix_;

  int rows_; 
  int cols_;
  
public:
  // creating empty matrix 
  DenseMatrix() : //! makes Matrix of zero length
    matrix_()
  { 
    rows_ = 0;
    cols_ = 0;
  }
  
  //! Copy Constructor
  DenseMatrix(const DenseMatrix<T> & org) : 
    matrix_(org.matrix_) , rows_(org.rows_) , cols_(org.cols_)
  {
  } 

  //! make matrix with 'rows' rows and 'cols' columns,
  //! maximum 'nz' non zero values in each row 
  //! and intialize all values with 'val'
  DenseMatrix(int rows, int cols)  
  { 
    resize(rows,cols);
  }

  void resize(int rows, int cols)
  { 
    matrix_.resize(rows);
    rows_ = rows; 
    cols_ = cols;
    
    for(int row=0; row<rows_; ++row) 
    {
      matrix_[row].resize(cols);
      for(register int col=0; col<cols_; ++col) matrix_[row][col] = 0;
    }
  }
  
/*******************************/
/*  Access and info functions  */
/*******************************/
  int rows() const {return rows_; }
  int cols() const {return cols_; }

  //int base() const {return base_;}
  T& operator() (int row, int col) 
  {
    assert( row>= 0); 
    assert( col>= 0);

    assert( row < rows_);
    assert( col < cols_);
    return matrix_[row][col];
  }        
  const T& operator() (int row, int col) const
  {
    assert( row >= 0); 
    assert( col >= 0);

    assert( row < rows_);
    assert( col < cols_);
    return matrix_[row][col];
  }        

  RowType & operator [] (int row) { return matrix_[row]; }
  const RowType & operator [] (int row) const { return matrix_[row]; }

  // result = this * vec 
  void mult(const T * vec, RowType & result) const 
  {
    assert( (int) result.size() == rows() );
    for(int row=0; row<rows(); ++row)
    {
      result[row] = 0;
      for(int col=0; col<cols(); ++col) 
      {
        result[row] += matrix_[row][col]*vec[col]; 
      }
    }
  }
    
  // result = this * vec 
  void multOEM(const T * vec, T * result) const 
  {
    assert( (int) result.size() == rows() );
    for(int row=0; row<rows(); ++row)
    {
      result[row] = 0;
      for(int col=0; col<cols(); ++col) 
      {
        result[row] += matrix_[row][col]*vec[col]; 
      }
    }
  }
    
  // result = this * vec 
  void mult(const RowType & vec, RowType & result) const 
  {
    this->mult(&vec[0],result);
  }
    
  // this = A * B
  void multiply(const DenseMatrix & A, const DenseMatrix & B)
  {
    assert( A.cols() == B.rows() );
   
    resize(A.rows() , B.cols());
    
    for(register int row=0; row<rows(); ++row)
    {
      for(register int col=0; col<cols(); ++col)
      {
        for(register int k=0; k<A.cols(); ++k)
          matrix_[row][col] += A[row][k] * B[k][col];
      }
    } 
  };

  DenseMatrix<T> & operator = (const DenseMatrix & org)
  {
    rows_ = org.rows_; 
    cols_ = org.cols_; 

    matrix_ = org.matrix_;
    return *this;
  } 
  
  DenseMatrix<T> & operator += (const DenseMatrix & org)
  {
    assert( rows() == rows() );
    assert( cols() == org.cols() );
    for(register int row=0; row<rows(); ++row)
    {
      for(register int col=0; col<cols(); ++col) 
        matrix_[row][col] += org.matrix_[row][col];
    }
    return *this;
  } 

  DenseMatrix<T> & operator -= (const DenseMatrix & org)
  {
    assert( rows() == rows() );
    assert( cols() == org.cols() );
    for(register int row=0; row<rows(); ++row)
    {
      for(register int col=0; col<cols(); ++col) 
        matrix_[row][col] -= org.matrix_[row][col];
    }
    return *this;
  } 

  void print(std::ostream & s) const
  {
    for(int i=0; i<rows(); ++i) 
    {
      for(int j=0; j<cols(); ++j) 
        s << matrix_[i][j] << " "; 
      s << std::endl;
    }
  }

  // set all matrix entries to zero 
  void clear() 
  {
    for(register int row=0; row<rows(); ++row)
    {
      for(register int col=0; col<cols(); ++col) 
        matrix_[row][col] = 0; 
    }
  }
};

//! Send vector to output stream
template<typename K>
std::ostream& operator<< (std::ostream& s, const DenseMatrix<K> & matrix)
{
  matrix.print(s);
  return s;
}

template <class T>
class BlockMatrix 
{
public: 
  typedef T Ttype;  //! remember the value type
  typedef BlockMatrix<T> ThisType;

private:
  typedef DenseMatrix<T> DenseMatrixType;
  typedef SparseRowMatrix< DenseMatrixType * > BlockMatrixType;

  typedef std::stack< DenseMatrixType * > MatrixStackType;

  DenseMatrixType * init_;
  BlockMatrixType matrix_;
  int localRows_;
  int localCols_;

  bool rowWise_; 

  mutable MatrixStackType freeStack_;
  
public:
  BlockMatrix(); //! makes Matrix of zero length
  BlockMatrix(const BlockMatrix<T> &S); //! Copy Constructor

  //! make matrix with 'rows' rows and 'cols' columns,
  //! maximum 'nz' non zero values in each row 
  //! and intialize all values with 'val'
  BlockMatrix(int rows, int cols, int lrows, int lcols, int nonZeros, bool
      rowWise = true) 
    : init_(0)
    , matrix_(rows,cols,nonZeros,init_,false)
    , localRows_(lrows),localCols_(lcols) 
    , rowWise_(rowWise)
  {}
  
  //! free memory for values_ and col_
  ~BlockMatrix()
  {
    clear();

    while( !freeStack_.empty() )
    {
      DenseMatrixType * dm = freeStack_.top();
      freeStack_.pop();
      delete dm; 
    }
  }
  
  DenseMatrixType & getMatrix() const 
  {
    if(freeStack_.empty())
    {
      return *(new DenseMatrixType(localRows_,localCols_)); 
    }
    else 
    {
      DenseMatrixType * dm = freeStack_.top();
      assert( dm->rows() == localRows_ );
      assert( dm->cols() == localCols_ );
      freeStack_.pop();
      return *dm; 
    }
  }

  void clear () 
  {
    const int s = matrix_.numberOfValues(); 
    for(register int i=0; i<s; ++i) 
    {
      DenseMatrixType * dm = matrix_.popValue(i);
      if(dm)
      {
        dm->clear();
        freeStack_.push(dm); 
      }
    }
  }
    
/*******************************/
/*  Access and info functions  */
/*******************************/
  int size(int i) const {return matrix_.size(i);};

  int littleRows() const { return localRows_;}
  int littleCols() const { return localCols_;}

  DenseMatrixType & operator() (int i, int j) 
  {
    DenseMatrixType * dm = matrix_(i,j);
    return *dm;
  }        
  
  const DenseMatrixType & operator() (int i, int j) const        
  {
    DenseMatrixType * dm = matrix_(i,j);
    return *dm;
  }        

  void set(int row, int col, DenseMatrixType & val)
  {
    int r = (rowWise_) ? row : col; 
    int c = (rowWise_) ? col : row; 
    
    assert( matrix_(r,c) == 0 ); 
    matrix_.set(r,c,&val);
  }
  
  void add(int row, int col, DenseMatrixType & val)
  {
    int r = (rowWise_) ? row : col; 
    int c = (rowWise_) ? col : row; 
    
    DenseMatrixType * dm = matrix_(r,c);
    if(!dm)
    {
      matrix_.set(r,c,&val);
      return ;
    }
    *(dm) += val;
  }

  void resize( int nsize) 
  {
    matrix_.clear();
    matrix_.resize(nsize);
  }
  
  void resize( int n, int m) 
  {
    matrix_.clear();
    matrix_.resize(n,m);
  }

#if 0
  void multiply(const ThisType & A , const ThisType & B )
  {
    assert( A.size(1) == B.size(0) );

    //matrix_.clear();
    
    std::cout << "Multiply: resize Matrix \n";
    resize(A.size(0), B.size(1));
    matrix_.clear();
    
    const int r = A.littleRows();
    const int c = B.littleCols();

    assert( A.matrix_.NumNonZeros() == B.matrix_.NumNonZeros() );

    DenseMatrixType tmp;

    for(int row=0; row<size(0); ++row)
    {
      for(int j=0; j<size(1); ++j)
      {
        DenseMatrixType dm(r,c);
        bool hasValue = false;
        for(int k=0; k<A.matrix_.NumNonZeros(); ++k)
        {
          std::pair < DenseMatrixType *, int> a = A.matrix_.realValue(row,k);
          if(!a.first) break;
          
          DenseMatrixType * b = B.matrix_ (a.second,j);
          if(!b) continue;

          tmp.multiply(*a.first,*b);
          dm += tmp;
          hasValue = true;
        }
        if(hasValue) set(row,j,dm);
      }
    }
    std::cout << "Done Multiply\n";
  }
#endif
  
  /*
  void multiply(const ThisType & A , const ThisType & B )
  {
    assert( A.size(1) == B.size(0) );

    //matrix_.clear();
    
    std::cout << "Multiply: resize Matrix \n";
    resize(A.size(0), B.size(1));
    matrix_.clear();
    
    const int r = A.littleRows();
    const int c = B.littleCols();

    DenseMatrixType tmp;

    for(int i=0; i<size(0); ++i)
    {
      for(int j=0; j<size(1); ++j)
      {
        DenseMatrixType dm(r,c);
        bool hasValue = false;
        for(int k=0; k<A.size(1); ++k)
        {
          DenseMatrixType * a = A.matrix_ (i,k);
          if(!a) continue;
          
          DenseMatrixType * b = B.matrix_ (k,j);
          if(!b) continue;

          tmp.multiply(*a,*b);
          dm += tmp;
          hasValue = true;
        }
        if(hasValue) set(i,j,dm);
      }
    }
    std::cout << "Done Multiply\n";
  }
  */

  // result = this * vec 
  void multOEM(const T * vec, T * result)
  {
    std::vector< T > v(localCols_); 
    std::vector< T > ret(localRows_); 

    for(register int r=0; r<size(0); ++r)
    {
      const int row = r * localRows_;

      // set right hand side to zero 
      for(register int k=0; k<localRows_; ++k) result[row+k] = 0;
      
      for(register int c=0; c<matrix_.NumNonZeros();++c)
      {
        std::pair< DenseMatrixType * , int > p = matrix_.realValue(r,c);
        if(p.first) 
        {
          const int col = p.second * localCols_;
          const T * v = &vec[col];
          p.first->mult(v,ret);

          for(register int k=0; k<localRows_; ++k)
            result[row+k] += ret[k];
        }
        else break;
      }
    }
  }
  
  // this += A 
  void add(const BlockMatrix<T> & A)
  {
    assert( matrix_.NumNonZeros() == A.matrix_.NumNonZeros() );
    assert( size(0) == A.size(0) );
    assert( size(1) == A.size(1) );
    
    const int nz = matrix_.NumNonZeros();
    const int vecSize = size(0) * nz; 
    
    for(register int r=0; r< vecSize; ++r)
    {
#ifndef NDEBUG
      std::pair< DenseMatrixType * , int > p = matrix_.realValue(r);
      std::pair< DenseMatrixType * , int > a = A.matrix_.realValue(r);
      // assert that column number is the same 
      assert( a.second == p.second );
#endif
      DenseMatrixType * m = matrix_.val(r);
      if(m) 
      {
        DenseMatrixType * a = A.val(r);
        assert( a );

        *m += *a; 
      }
    }
  }

  void print(std::ostream & s) const 
  {
    for(int i=0; i<matrix_.size(0); ++i) 
    {
      for(int j=0; j<matrix_.size(1); ++j) 
      {
        DenseMatrixType * m = matrix_(i,j);
        if(m) 
        {
          s << (*m); 
        }
      }
      s << std::endl;
    }
  }

};

//! Send vector to output stream
template<typename K>
std::ostream& operator<< (std::ostream& s, const BlockMatrix<K> & matrix)
{
  matrix.print(s);
  return s;
}

} // end namespace Dune 
#endif  
