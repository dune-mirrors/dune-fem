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
    if( (rows == rows_) && (cols == cols_) ) return ;
    
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
    assert( ((int) result.size() != rows()) ? 
        (std::cout << result.size() << " s|r " << rows() << "\n",0) : 1 );
    const int nRow= rows();
    const int nCol= cols();
    for(int row=0; row<nRow; ++row)
    {
      result[row] = 0;
      for(int col=0; col<nCol; ++col) 
      {
        result[row] += matrix_[row][col]*vec[col]; 
      }
    }
  }
    
  // result = this * vec 
  void multOEM(const T * vec, T * result) const 
  {
    const int nRow= rows();
    const int nCol= cols();
    for(int row=0; row<nRow; ++row)
    {
      result[row] = 0;
      for(int col=0; col<nCol; ++col) 
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
    
  // result = this^T * vec 
  void multTransposed(const RowType & vec, RowType & result) const 
  {
    assert( (int) result.size() == cols() );
    const int nCols = cols();
    const int nRows = rows();
    for(int col=0; col<nCols; ++col) 
    {
      result[col] = 0;
      for(int row=0; row<nRows; ++row)
      {
        result[col] += matrix_[row][col]*vec[row]; 
      }
    }
  }
    
  // this = A * B
  void multiply(const DenseMatrix & A, const DenseMatrix & B)
  {
    assert( A.cols() == B.rows() );
   
    resize( A.rows() , B.cols() );
    
    const int nRows = rows();
    const int nCols = cols();
    const int Acols = A.cols();
    for(int row=0; row<nRows; ++row)
    {
      for(int col=0; col<nCols; ++col)
      {
        T sum = 0;
        for(int k=0; k<Acols; ++k)
        {
          sum += A[row][k] * B[k][col];
        }
        matrix_[row][col] = sum; 
      }
    } 
  };

  // this = A * B
  void multiplyTransposed(const DenseMatrix & A, const DenseMatrix & B)
  {
    assert( A.cols() == B.cols() );
   
    resize(A.rows() , B.rows());
    
    for(int row=0; row<rows(); ++row)
    {
      for(int col=0; col<cols(); ++col)
      {
        T sum = 0;
        for(int k=0; k<A.cols(); ++k)
        {
          sum += A[row][k] * B[col][k];
        }
        matrix_[row][col] = sum; 
      }
    } 
  };

  // this = A * B
  void multiply_AT_A(const DenseMatrix & A) 
  {
    resize(A.cols() , A.cols());
    
    for(int row=0; row<rows(); ++row)
    {
      for(int col=0; col<cols(); ++col)
      {
        T sum = 0;
        const int aRows = A.rows();
        for(int k=0; k<aRows; ++k)
        {
          sum += A[k][row] * A[k][col];
        }
        matrix_[row][col] = sum; 
      }
    } 
  }

  DenseMatrix<T> & operator = (const DenseMatrix & org)
  {
    rows_ = org.rows_; 
    cols_ = org.cols_; 

    matrix_ = org.matrix_;
    return *this;
  } 
  
  DenseMatrix<T> & operator += (const DenseMatrix & org)
  {
    const int nRows = rows();
    const int nCols = cols();
    assert( nRows == org.rows() );
    assert( nCols == org.cols() );
    for(int row=0; row<nRows; ++row)
    {
      for(int col=0; col<nCols; ++col) 
      {
        matrix_[row][col] += org.matrix_[row][col];
      }
    }
    return *this;
  } 

  DenseMatrix<T> & operator -= (const DenseMatrix & org)
  {
    assert( rows() == rows() );
    assert( cols() == org.cols() );
    for(int row=0; row<rows(); ++row)
    {
      for(int col=0; col<cols(); ++col) 
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
    for(int row=0; row<rows(); ++row)
    {
      for(int col=0; col<cols(); ++col) 
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

  BlockMatrixType matrix_;
  int localRows_;
  int localCols_;

  bool rowWise_; 

  mutable MatrixStackType freeStack_;

  int nonZeros_;

public:
  //! makes Matrix of zero length
  BlockMatrix() : matrix_(), localRows_(0), localCols_(0) ,
    rowWise_(true), nonZeros_(0) {} 
  BlockMatrix(const BlockMatrix<T> &S); //! Copy Constructor

  //! make matrix with 'rows' rows and 'cols' columns,
  //! maximum 'nz' non zero values in each row 
  //! and intialize all values with 'val'
  BlockMatrix(int rows, int cols, int lrows, 
              int lcols, int nonZeros, bool rowWise = true) 
  {
    reserve(rows,cols,lrows,lcols,nonZeros,rowWise);
  }
  
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

  void reserve(int rows, int cols, int lrows, int lcols, 
               int nonZeros, bool rowWise = true) 
  {
    DenseMatrixType * init = 0;

    rowWise_ = rowWise;

    localRows_ = (rowWise_) ? lrows : lcols; 
    localCols_ = (rowWise_) ? lcols : lrows; 
    nonZeros_ = nonZeros;

    int r = (rowWise_) ? rows : cols; 
    int c = (rowWise_) ? cols : rows; 
    matrix_.reserve(r,c,nonZeros,init);
  }
  
  void freeMatrix(DenseMatrixType & dm) const 
  {
    dm.clear();
    freeStack_.push(&dm); 
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
        freeMatrix(*dm);
      }
    }
  }
    
/*******************************/
/*  Access and info functions  */
/*******************************/
  int size(int i) const {return matrix_.size(i);};

  int littleRows() const { return localRows_;}
  int littleCols() const { return localCols_;}

  DenseMatrixType & operator() (int row, int col) 
  {
    int r = (rowWise_) ? row : col; 
    int c = (rowWise_) ? col : row; 
    
    DenseMatrixType * dm = matrix_(r,c);
    return *dm;
  }        
  
  const DenseMatrixType & operator() (int row, int col) const        
  {
    int r = (rowWise_) ? row : col; 
    int c = (rowWise_) ? col : row; 
    
    DenseMatrixType * dm = matrix_(r,c);
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
    resize( nsize, nsize );
  }

  void resize( int rows, int cols ) 
  {
    matrix_.resize( rows , cols ); 
  }

  void resort () 
  {
    const int rows = matrix_.size(0); 
    for (int i =0 ; i<rows; ++i) 
    {
      resortRow(i);
    }
  }

  void resortRow ( int row ) 
  {
    matrix_.resortRow( row );
  }
  
  void multiply(const ThisType & A , const ThisType & B )
  {
    assert( A.size(1) == B.size(0) );
    // A has to be rowWise oriented matrix 
    assert( A.rowWise_ == true );
    // B has to be colWise oriented matrix 
    //assert( B.rowWise_ == false );

    std::cout << "Multiply: resize Matrix \n";
    clear();

    //if( nonZeros_ == 0) nonZeros_ = 2*std::max(A.nonZeros_, B.nonZeros_ ); 
    if( nonZeros_ == 0) nonZeros_ = std::max(A.nonZeros_, B.nonZeros_ ); 

    // reserve memory and set right new block sizes
    reserve( A.size(0), B.size(1), 
             A.littleRows(), 
             B.littleCols(), 
             nonZeros_  );
    
    DenseMatrixType tmp;
    for(int row=0; row<size(0); ++row)
    {
      for(int j=0; j<size(1); ++j)
      {
        // get empty matrix from stack 
        DenseMatrixType & dm = getMatrix();

        bool hasValue = false;
        // get number of nonZeros of row 
        const int nonZeros = A.matrix_.NumNonZeros(row);
        for(int k=0; k<nonZeros; ++k)
        {
          std::pair < DenseMatrixType *, int> a = A.matrix_.realValue(row,k);
          // use A.resort to achieve this feature 
          assert( a.first );
          //if(!a.first) break;
          
          // serach corresponding matrix entry in B  
          DenseMatrixType * b = B.matrix_ (a.second,j);
          // if not available, continue 
          if(!b) continue;
          
          tmp.multiply(*a.first,*b);
          dm += tmp;
          hasValue = true;
        }
        
        // if value found, set entry in matrix 
        if(hasValue) 
        {
          set(row,j,dm);
        }
        // else free memory 
        else 
        {
          // clears dm and put to stack 
          freeMatrix(dm);
        }
      }
    }
    //print(std::cout);
    std::cout << "Done Multiply\n";
  }

  void add(const ThisType & A)
  {
    std::cout << "Start add of Matrices \n";
    assert( A.size(0)  == this->size(0) );
    assert( A.size(1)  == this->size(1) );
    assert( A.rowWise_ == rowWise_);

    const int nRows = size(0);
    for(int row=0; row<nRows; ++row)
    {
      const int nonZeros = A.matrix_.NumNonZeros(row); 
      for(int k=0; k<nonZeros; ++k)
      {
        std::pair < DenseMatrixType *, int> a = A.matrix_.realValue(row,k);
        assert( a.first );

        // serach corresponding matrix entry in B  
        DenseMatrixType * b = matrix_(row,a.second);
        // if not available, continue 
        if(!b)
        {
          // get empty matrix from stack 
          DenseMatrixType & dm = getMatrix();
          dm = *a.first;
          set(row,a.second,dm);
        }
        else 
        {
          *b += (*a.first);
        }
      }
    }
    std::cout << "Done Add\n";
  }

  // result = this * vec 
  void multOEM(const T * vec, T * result) const
  {
    std::vector< T > v(localCols_); 
    std::vector< T > ret(localRows_); 

    for(int r=0; r<size(0); ++r)
    {
      const int row = r * localRows_;

      // set right hand side to zero 
      for(int k=0; k<localRows_; ++k) result[row+k] = 0;
      
      const int nonZeros = matrix_.NumNonZeros(r);
      for(int c=0; c<nonZeros; ++c)
      {
        std::pair< DenseMatrixType * , int > p = matrix_.realValue(r,c);
        assert( p.first );
        const int col = p.second * localCols_;
        const T * v = &vec[col];
        p.first->mult(v,ret);

        for(register int k=0; k<localRows_; ++k)
          result[row+k] += ret[k];
      }
    }
  }
  
  // result = this * vec 
  void multOEMAdd(const T * vec, T * result) const
  {
    std::vector< T > v(localCols_); 
    std::vector< T > ret(localRows_); 

    for(register int r=0; r<size(0); ++r)
    {
      const int row = r * localRows_;

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
  
  //! this precondition is from right side 
  bool rightPrecondition () const { return true; }
  
  // result = this * vec 
  void precondition(const T * u, T * x) const
  {
    const int nRows = matrix_.size(0);
    const double omega = 1.1;

    const T * uox = x;
    
    // (D - omega E) x = x_old (=u)  
    for(int r=0; r<nRows; ++r)
    {
      const int nonZeros = matrix_.NumNonZeros(r);
      int row = r * localRows_;
      for(int i=0; i<localRows_; ++i)
      {
        double dot = 0.0;
        double diag = 1.0;

        for(int k=0; k<nonZeros; ++k)
        {
          std::pair < DenseMatrixType *, int> a = matrix_.realValue(r,k);
          assert( a.first );

          DenseMatrixType & dm = *a.first;
          if( a.second < r ) 
          {
            int col = a.second * localCols_;
            for(int c=0; c<localCols_; ++c)
            {
              dot += dm[i][c] * uox[col];
              ++col;
            }
          }
          if( a.second == r ) 
          {
            int col = a.second * localCols_;
            for(int c=0; c<i; ++c)
            {
              dot += dm[i][c] * uox[col];
              ++col;
            }
            diag = dm[i][i];
          }
        }
          
        x[row] = (u[row] - omega*dot) / diag;
        ++row;
      }
    }

    // (D - omega E) x = x_old (=u)  
    for(int r=nRows-1; r>= 0; --r)
    {
      const int nonZeros = matrix_.NumNonZeros(r);
      for(int i=localRows_-1; i>=0; --i)
      {
        int row = (r * localRows_) +i;
        double dot = 0.0;
        double diag = 1.0;

        for(int k=0; k<nonZeros; ++k)
        {
          std::pair < DenseMatrixType *, int> a = matrix_.realValue(r,k);
          DenseMatrixType & dm = *a.first;
          if( a.second > r ) 
          {
            int col = a.second * localCols_;
            for(int c=0; c<localCols_; ++c)
            {
              dot += dm[i][c] * uox[col];
              ++col;
            }
          }
          if( a.second == r ) 
          {
            int col = (a.second * localCols_) + i+1;
            for(int c=i+1; c<localCols_; ++c)
            {
              dot += dm[i][c] * uox[col];
              ++col;
            }
            diag = dm[i][i];
          }
        }
         
        x[row] -= omega * dot / diag; 
      }
    }
  }
  
  // result = this * vec 
  void setDiag(const T * diag)
  {
    const int nRows = size(0);
    for(int r=0; r<nRows; ++r)
    {
      const int row = r * localRows_;

      const int nonZeros = matrix_.NumNonZeros(r); 
      for(int c=0; c<nonZeros; ++c)
      {
        std::pair< DenseMatrixType * , int > p = matrix_.realValue(r,c);
        assert( p.first );
        const int littleR = p.first->rows();
        DenseMatrixType & dm = *p.first;
        for(int i=0; i<littleR; ++i)
        {
          double val = diag[row+i];
          const int littleC = dm.cols();
          for(int lc = 0; lc<littleC; ++lc)
          {
            dm[i][lc] *= val; 
          }
        }
      }
    }
  }

  void getDiag(T * diag) 
  {
    const int nRows = size(0); 
    for(int r=0; r<nRows; ++r)
    {
      const int row = r * localRows_;
      
      DenseMatrixType * d = matrix_(r,r);
      assert( d );
      DenseMatrixType & dm = *d; 
      
      // set right hand side to zero 
      for(int k=0; k<localRows_; ++k) diag[row+k] = dm[k][k];
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

template <class RowSpaceType, class ColumnSpaceType> 
class BlockMatrixObject
{
  typedef typename RowSpaceType::GridType::template Codim<0>::Entity EntityType;

  typedef BlockMatrixObject<RowSpaceType,ColumnSpaceType> ThisType;
public:  

  typedef BlockMatrix<double> MatrixType;
  typedef MatrixType PreconditionMatrixType;
    
  template <class MatrixImp> 
  class LocalMatrix
  {
    typedef MatrixImp MatrixType;
    typedef DenseMatrix<double> LittleBlockType;

    MatrixType & matrix_; 
    
    const int rowIndex_;
    const int colIndex_;

    LittleBlockType & localMatrix_;

  public:  
    LocalMatrix(MatrixType & m,
                const EntityType& rowEntity,
                const RowSpaceType & rowSpace,
                const EntityType& colEntity,
                const ColumnSpaceType & colSpace)
      : matrix_(m)
      , rowIndex_(rowSpace.indexSet().index(rowEntity))
      , colIndex_(colSpace.indexSet().index(colEntity))
      , localMatrix_(matrix_.getMatrix())
    {
    }

    // add block to block matrix 
    ~LocalMatrix()
    {
      // finalize by adding block to global matrix 
      addLocalMatrix();
    }

  private: 
    LocalMatrix(const LocalMatrix &);

    void addLocalMatrix() 
    {
      if( matrix_.size(0) > 0 )
      {
        matrix_.add(rowIndex_,colIndex_, localMatrix_ );
      }
    }

  public:
    int rows () const { return localMatrix_.rows(); }
    int cols () const { return localMatrix_.cols(); }

    void add(int localRow, int localCol , const double value)
    {
      assert( localRow >= 0 );
      assert( localCol >= 0 );

      assert( localRow < localMatrix_.rows() );
      assert( localCol < localMatrix_.cols() );

      localMatrix_[localRow][localCol] += value; 
    }

    double get(int localRow, int localCol ) const 
    {
      assert( localRow >= 0 );
      assert( localCol >= 0 );

      assert( localRow < localMatrix_.rows() );
      assert( localCol < localMatrix_.cols() );

      return localMatrix_[localRow][localCol]; 
    }

    //! set matrix enrty to value 
    void set(int localRow, int localCol, const double value)
    {
      assert( localRow >= 0 );
      assert( localCol >= 0 );

      assert( localRow < localMatrix_.rows() );
      assert( localCol < localMatrix_.cols() );

      localMatrix_[localRow][localCol] = value; 
    }

    //! clear all entries belonging to local matrix 
    void clear ()
    {
      localMatrix_.clear();
    }

    //! resort all global rows of matrix to have ascending numbering 
    void resort ()
    {
      //assert( matrix_.size(0) > 0 );
      //matrix_.resortRow( rowIndex_ );
    }
  };

  
public:
  typedef LocalMatrix<MatrixType> LocalMatrixType;

  const RowSpaceType & rowSpace_; 
  const ColumnSpaceType & colSpace_;

  int size_;
  
  int numRowBaseFct_;
  int numColBaseFct_;

  MatrixType matrix_; 
  const bool preconditioning_;

  //! setup matrix handler 
  BlockMatrixObject(const RowSpaceType & rowSpace, 
                    const ColumnSpaceType & colSpace,
                    bool preconditioning) 
    : rowSpace_(rowSpace)
    , colSpace_(colSpace) 
    , size_(-1)
    , numRowBaseFct_(-1)
    , numColBaseFct_(-1)
    , matrix_()
    , preconditioning_(preconditioning)
  {
    assert( rowSpace_.indexSet().size(0) == 
            colSpace_.indexSet().size(0) ); 
    reserve(true);
  }

  //! return reference to stability matrix 
  MatrixType & matrix() { return matrix_; }

  //! resize all matrices and clear them 
  void clear() 
  {
    matrix_.clear();
  }

  //! return true if precoditioning matrix is provided 
  bool hasPcMatrix () const { return preconditioning_; }
  PreconditionMatrixType& pcMatrix () { return matrix_; }

  //! resize all matrices and clear them 
  void resize(bool verbose = false) 
  {
    if( ! hasBeenSetup() ) 
    {
      reserve(); 
    }
    else 
    {
      size_ = rowSpace_.indexSet().size(0); 
      if(verbose)
      {
        std::cout << "Resize Matrix with (" << size_ << "," << size_ << ")\n";
      }
      matrix_.resize(size_);
    }
  }

  //! returns true if memory has been reserved
  bool hasBeenSetup () const 
  {
    return (numRowBaseFct_ > 0);
  }

  //! reserve memory corresponnding to size of spaces 
  void reserve(bool verbose = false ) 
  {
    // if empty grid do nothing (can appear in parallel runs)
    if( (rowSpace_.begin() != rowSpace_.end()) && 
        (colSpace_.begin() != colSpace_.end()) )
    {
      // get number of elements 
      size_ = rowSpace_.indexSet().size(0); 

      numRowBaseFct_ = rowSpace_.getBaseFunctionSet(*(rowSpace_.begin())).numBaseFunctions();
      numColBaseFct_ = colSpace_.getBaseFunctionSet(*(colSpace_.begin())).numBaseFunctions();

      if(verbose) 
      {
        std::cout << "Reserve Matrix with (" << size_ << "," << size_ << ")\n";
        std::cout << "Number of base functions = (" << numRowBaseFct_ << "," << numColBaseFct_ << ")\n";
      }

      // factor for non-conforming grid is 4 in 3d and 2 in 2d  
      //const int factor = (Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));

      // upper estimate for number of neighbors 
      enum { dim = RowSpaceType :: GridType :: dimension };

      // number of neighbors + 1 
      const int factor = (dim * 2) + 1; 

      // upper estimate for number of neighbors 
      //enum { dim = RowSpaceType :: GridType :: dimension };
      //rowMaxNumbers_ *= (factor * dim * 2) + 1; // e.g. 7 for dim = 3
      matrix_.reserve(size_,size_,numRowBaseFct_,numColBaseFct_,factor);
    }
  }

  //! mult method of matrix object used by oem solver
  void multOEM(const double * arg, double * dest) const 
  {
    matrix_.multOEM(arg,dest);
  }

  //! resort row numbering in matrix to have ascending numbering 
  void resort() 
  {
    matrix_.resort();
  }

  //! empty method as we use right preconditioning here
  void createPreconditionMatrix() {}
};




} // end namespace Dune 
#endif  
