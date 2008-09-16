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
  typedef std::vector < T > RowType;

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

  //! this = A^T * A 
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

  //! set all values of the matrix to given value 
  DenseMatrix<T> & operator = (const T & val)
  {
    for(int row=0; row<rows(); ++row)
    {
      for(int col=0; col<cols(); ++col)
      {
        matrix_[row][col] = val;
      }
    }
    return *this;
  } 
  
  //! set all values of the matrix to values of given matrix 
  DenseMatrix<T> & operator = (const DenseMatrix & org)
  {
    rows_ = org.rows_; 
    cols_ = org.cols_; 

    matrix_ = org.matrix_;
    return *this;
  } 
  
  //! add matrix 
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

  //! substract matrix 
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

  //! print matrix 
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
    // clear all matrices to stack 
    clear (lrows, lcols);

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

  void clear (const int lrows, const int lcols) 
  {
    if( lrows != localRows_ || lcols != localCols_ )
    {
      const int s = matrix_.numberOfValues(); 
      for(register int i=0; i<s; ++i) 
      {
        DenseMatrixType * dm = matrix_.popValue(i);
        if(dm)
        {
          dm->resize(lrows,lcols);
          freeMatrix(*dm);
        }
      }
    }
    else 
    {
      clear();
    }
  }
    
/*******************************/
/*  Access and info functions  */
/*******************************/
  int size(int i) const {return matrix_.size(i);};

  int littleRows() const { return localRows_;}
  int littleCols() const { return localCols_;}

  DenseMatrixType& createEntry(int row, int col)
  {
    DenseMatrixType& dm = getMatrix();
    set(row,col,dm); 
    return dm;
  }
  
  DenseMatrixType & operator() (int row, int col) 
  {
    int r = (rowWise_) ? row : col; 
    int c = (rowWise_) ? col : row; 
    
    DenseMatrixType * dm = matrix_(r,c);
    return (dm) ? *dm : createEntry(r,c);
  }        

  const DenseMatrixType & operator() (int row, int col) const        
  {
    int r = (rowWise_) ? row : col; 
    int c = (rowWise_) ? col : row; 
    
    DenseMatrixType * dm = matrix_(r,c);
    return (dm) ? *dm : createEntry(r,c);
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
        const int nonZeros = A.matrix_.numNonZeros(row);
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
      const int nonZeros = A.matrix_.numNonZeros(row); 
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
      
      const int nonZeros = matrix_.numNonZeros(r);
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

      for(register int c=0; c<matrix_.numNonZeros();++c)
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
      const int nonZeros = matrix_.numNonZeros(r);
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
      const int nonZeros = matrix_.numNonZeros(r);
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

      const int nonZeros = matrix_.numNonZeros(r); 
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

template <class DomainSpace, class RangeSpace, class OperatorTraits >
class BlockMatrixObject;

template <class DomainSpaceImp, class ColSpaceImp = DomainSpaceImp>
struct BlockMatrixTraits
{
  typedef DomainSpaceImp RowSpaceType;
  typedef ColSpaceImp ColumnSpaceType;
  typedef BlockMatrixTraits<RowSpaceType,ColumnSpaceType> ThisType;

  template <class OperatorTraits>
  struct MatrixObject
  {
    typedef BlockMatrixObject<RowSpaceType,ColumnSpaceType,OperatorTraits> MatrixObjectType;
  };
};


//! matrix object holding a blockamtrix
template <class DomainSpaceImp, class RangeSpaceImp, class OperatorTraits> 
class BlockMatrixObject
{
public:  
  typedef DomainSpaceImp DomainSpaceType;
  typedef RangeSpaceImp RangeSpaceType ;

  typedef typename OperatorTraits :: StencilType StencilType;


  //! number of rows of blocks 
  enum { littleRows = DomainSpaceType :: localBlockSize };
  //! number of columns of blocks 
  enum { littleCols = RangeSpaceType :: localBlockSize };

  typedef typename DomainSpaceType::GridType::template Codim<0>::Entity EntityType;

  typedef BlockMatrixObject<DomainSpaceType,RangeSpaceType,OperatorTraits> ThisType;

  typedef BlockMatrix<double> MatrixType;
  typedef MatrixType PreconditionMatrixType;

  template <class MatrixObjectImp> 
  class LocalMatrix;

  struct LocalMatrixTraits
  {
    typedef DomainSpaceImp DomainSpaceType ;
    typedef RangeSpaceImp RangeSpaceType;
    typedef typename DomainSpaceImp :: RangeFieldType RangeFieldType;
    typedef LocalMatrix<ThisType> LocalMatrixType;
    typedef DenseMatrix<RangeFieldType> LittleBlockType;
  };

  template <class MatrixObjectImp> 
  class LocalMatrix : public LocalMatrixDefault<LocalMatrixTraits>
  {
  public:  
    //! type of base class 
    typedef LocalMatrixDefault<LocalMatrixTraits> BaseType;

    //! type of matrix object 
    typedef MatrixObjectImp MatrixObjectType;
    //! type of matrix 
    typedef typename MatrixObjectImp :: MatrixType MatrixType;
    //! type of entries of little blocks 
    typedef typename DomainSpaceType :: RangeFieldType DofType;

    //! type of little blocks 
    typedef DenseMatrix<DofType> LittleBlockType;

    //! type of row mapper 
    typedef typename MatrixObjectType :: RowMapperType RowMapperType;
    //! type of col mapper 
    typedef typename MatrixObjectType :: ColMapperType ColMapperType;

  private:  
    const MatrixObjectType & matrixObj_; 
    const RowMapperType& rowMapper_;
    const ColMapperType& colMapper_;
    
    std::vector< int > rows_;
    std::vector< int > cols_;
    std::vector< std::vector<LittleBlockType*> > matrices_;

  public:  
    LocalMatrix(const MatrixObjectType & mObj,
                const DomainSpaceType & rowSpace,
                const RangeSpaceType & colSpace,
                const RowMapperType & rowMapper,
                const ColMapperType & colMapper)
      : BaseType(rowSpace,colSpace) 
      , matrixObj_(mObj)
      , rowMapper_(rowMapper)
      , colMapper_(colMapper)
      , matrices_() 
    {
    }

    //! initialize local matrix to entities 
    void init(const EntityType& rowEntity,
              const EntityType& colEntity)
    {
      // initialize base functions sets 
      BaseType :: init ( rowEntity , colEntity );

      MatrixType& matrix = matrixObj_.matrix();
      // get global block numbers 
      const size_t rows = rowMapper_.numEntityDofs(rowEntity);
      {
        rows_.resize( rows );
        for(size_t i=0; i<rows; ++i)
        {
          rows_[i] = rowMapper_.mapToGlobal( rowEntity, i );
        }
      }

      // get global block numbers 
      const size_t cols = colMapper_.numEntityDofs(colEntity);
      {
        cols_.resize( cols );
        for(size_t i=0; i<cols; ++i)
        {
          cols_[i] = colMapper_.mapToGlobal( colEntity, i );
        }
      }

      matrices_.resize( rows );
      for(size_t i=0; i<rows; ++i)
      {
        matrices_[i].resize( cols );
        for(size_t j=0; j<cols; ++j) 
        {
          // get pointer to block matrices 
          matrices_[i][j] = &matrix(rows_[i],cols_[j]);
        }
      }
    }

  private: 
    //! prohibit copy cosntructor 
    LocalMatrix(const LocalMatrix &);

    // check whether given (row,col) pair is valid
    void check(int localRow, int localCol) const
    {
      const size_t row = (int) localRow / littleRows;
      const size_t col = (int) localCol / littleCols;
      const int lRow = localRow%littleRows;
      const int lCol = localCol%littleCols;
      assert( row < matrices_.size() ) ;
      assert( col < matrices_[row].size() );
      assert( lRow < littleRows );
      assert( lCol < littleCols );
    }
    
    // return reference to entry 
    DofType& getValue(const int localRow, const int localCol) 
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      const int row = (int) localRow / littleRows;
      const int col = (int) localCol / littleCols;
      const int lRow = localRow%littleRows;
      const int lCol = localCol%littleCols;
      assert( matrices_[row][col] );
      return (*matrices_[row][col])[lRow][lCol];
    }

    // return reference to entry 
    const DofType& getValue(const int localRow, const int localCol) const 
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      const int row = (int) localRow / littleRows;
      const int col = (int) localCol / littleCols;
      const int lRow = localRow%littleRows;
      const int lCol = localCol%littleCols;
      assert( matrices_[row][col] );
      return (*matrices_[row][col])[lRow][lCol];
    }
  public:
    //! add value to matrix 
    void add(int localRow, int localCol , const DofType value)
    {
      getValue(localRow,localCol) += value; 
    }

    //! return matrix entry 
    const DofType get(int localRow, int localCol ) const 
    {
      return getValue(localRow,localCol); 
    }

    //! set matrix enrty to value 
    void set(int localRow, int localCol, const DofType value)
    {
      getValue(localRow,localCol) = value; 
    }

    //! set matrix enrty to value 
    void unitRow(const int localRow) 
    {
      const int col = this->columns();
      for(int localCol=0; localCol<col; ++localCol)
      {
        getValue(localRow,localCol) = 0;
      }
      getValue(localRow,localRow) = 1;
    }

    //! clear all entries belonging to local matrix 
    void clear ()
    {
      const size_t rows = rows_.size();
      const size_t cols = cols_.size();
      for(size_t i=0; i<rows; ++i) 
      {
        for(size_t j=0; j<cols; ++j) 
        {
          matrices_[i][j]->clear(); 
        }
      }
    }

    //! resort all global rows of matrix to have ascending numbering 
    void resort ()
    {
      MatrixType& matrix = matrixObj_.matrix();
      assert( matrix.size(0) > 0 );
      const size_t rows = rows_.size();
      for(size_t i=0; i<rows; ++i) 
      {
        matrix.resortRow( rows_[i] );
      }
    }
  };

public:
  //! type of local matrix 
  typedef LocalMatrix<ThisType> ObjectType;
  typedef ThisType LocalMatrixFactoryType;
  typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
  //! type of local matrix 
  typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

  typedef typename DomainSpaceType :: BlockMapperType RowMapperType; 
  typedef typename RangeSpaceType :: BlockMapperType ColMapperType; 

  typedef AdaptiveDiscreteFunction<DomainSpaceType> DestinationType;

  // row space 
  const DomainSpaceType & domainSpace_; 
  // column space 
  const RangeSpaceType & rangeSpace_;

  // sepcial row mapper 
  RowMapperType& rowMapper_;
  // special col mapper 
  ColMapperType& colMapper_;

  // number of grid sequence 
  int sequence_;
  
  // matrix 
  mutable MatrixType* matrix_; 
  // preconditioning flag 
  bool preconditioning_;

  // local matrix stack 
  mutable LocalMatrixStackType localMatrixStack_;

  //! constructor 
  //! \param rowSpace space defining row structure 
  //! \param colSpace space defining column structure 
  //! \param paramfile parameter file to read variables 
  //!         - Preconditioning: {0 == no, 1 == yes} used is SSOR 
  BlockMatrixObject(const DomainSpaceType & rowSpace, 
                    const RangeSpaceType & colSpace,
                    const std::string paramfile = "" ) 
    : domainSpace_(rowSpace)
    , rangeSpace_(colSpace) 
    , rowMapper_( rowSpace.blockMapper() )
    , colMapper_( colSpace.blockMapper() )
    , sequence_(-1)
    , matrix_(0)
    , preconditioning_(false)
    , localMatrixStack_(*this)
  {
    if( paramfile != "" )
    {
      int precon = 0;
      readParameter(paramfile,"Preconditioning",precon);
      preconditioning_ = (precon > 0) ? true : false;
    } 
    assert( domainSpace_.indexSet().size(0) == 
            rangeSpace_.indexSet().size(0) ); 
  }

  //! destructor 
  ~BlockMatrixObject() 
  {
    delete matrix_;
  }

  //! return reference to stability matrix 
  MatrixType & matrix() const 
  { 
    assert( matrix_ );
    return *matrix_; 
  }

  //! set all matrix entries to zero  
  void clear() 
  {
    matrix().clear();
  }

  //! return true if precoditioning matrix is provided 
  bool hasPreconditionMatrix () const { return preconditioning_; }

  //! return reference to preconditioner (here also systemMatrix)
  const PreconditionMatrixType& preconditionMatrix () const { return matrix(); }

  //! reserve memory corresponnding to size of spaces 
  void reserve(bool verbose = false ) 
  {
    if(sequence_ != domainSpace_.sequence() )
    {
      if( matrix_ ) 
      {
        delete matrix_;
        matrix_ = 0;
      }

      if( !matrix_ )
      {
        matrix_ = new MatrixType ();
      }

      // if empty grid do nothing (can appear in parallel runs)
      if( (domainSpace_.begin() != domainSpace_.end()) && 
          (rangeSpace_.begin() != rangeSpace_.end()) )
      {
        // get number of elements 
        int rowSize = rowMapper_.size();
        int colSize = colMapper_.size();

        // if verbose mode output information 
        if(verbose) 
        {
          std::cout << "Reserve Matrix with (" << rowSize << "," << colSize << ")\n";
          std::cout << "Number of base functions = (" << littleRows << "," << littleCols << ")\n";
        }

        // get number of non-zeros 
        const int nonZerosBlocks = 
          StencilType :: nonZerosEstimate( rangeSpace_ ) / RangeSpaceType :: localBlockSize ; 

        // upper estimate for number of neighbors 
        matrix().reserve(rowSize, colSize , littleRows , littleCols , nonZerosBlocks );
      }
      // set new sequence number 
      sequence_ = domainSpace_.sequence();
    }
  }

  //! apply matrix to discrete function
  template< class DomainFunction, class RangeFunction >
  void apply ( const DomainFunction &arg, RangeFunction &dest ) const
  {
    // apply matrix vector multiplication 
    matrix().multOEM( arg.leakPointer(), dest.leakPointer() );
    
    // communicate data 
    rangeSpace_.communicate( dest );
  }

  //! mult method of matrix object used by oem solver
  void multOEM(const double * arg, double * dest) const 
  {
    typedef AdaptiveDiscreteFunction< DomainSpaceType > DomainFunctionType;
    typedef AdaptiveDiscreteFunction< RangeSpaceType > RangeFunctionType;

    DomainFunctionType farg( "multOEM arg", domainSpace_, arg );
    RangeFunctionType fdest( "multOEM dest", rangeSpace_, dest );

    // matrix multiplication 
    apply(farg,fdest);
  }

  //! resort row numbering in matrix to have ascending numbering 
  void resort() 
  {
    matrix().resort();
  }

  //! print matrix 
  void print(std::ostream & s) const 
  {
    matrix().print(s);
  }

 //! interface method from LocalMatrixFactory 
  ObjectType* newObject() const
  {
    return new ObjectType(*this,
                          domainSpace_,
                          rangeSpace_,
                          rowMapper_,
                          colMapper_);
  }

  //! return local matrix 
  LocalMatrixType localMatrix(const EntityType& rowEntity,
                              const EntityType& colEntity) const
  {
    return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
  }

  //! empty method as we use right preconditioning here
  void createPreconditionMatrix() {}
};




} // end namespace Dune 
#endif  
