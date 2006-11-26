#ifndef DUNE_SPMATRIX_HH
#define DUNE_SPMATRIX_HH

#include <vector>

namespace Dune
{

//*****************************************************************
//
//  --SparseRowMatrix
//  
//! Compressed row sparse matrix, where only the nonzeros of a row a 
//! keeped
//*****************************************************************
template <class T>
class SparseRowMatrix 
{
public: 
  typedef T Ttype;  //! remember the value type

  typedef SparseRowMatrix<T> ThisType;

  enum { defaultCol = -1 };
  enum { firstCol = defaultCol + 1 };

private:
  T* values_ ;      //! data values (nz_ elements)
  int* col_;        //! row_ptr (dim_[0]+1 elements)
  int* nonZeros_;   //! row_ptr (dim_[0]+1 elements)
  int dim_[2];      //! dim_[0] x dim_[1] Matrix
  int nz_;          //! number of nonzeros per row

  int memSize_;
  bool sorted_;

  // temporary mem for resort 
  std::vector<int> newIndices_;
  std::vector<T> newValues_;
  
public:
  //! makes Matrix of zero length
  SparseRowMatrix(); 
  //! Copy Constructor
  SparseRowMatrix(const SparseRowMatrix<T> &S); 

  //! make matrix with 'rows' rows and 'cols' columns,
  //! maximum 'nz' non zero values in each row 
  //! and intialize all values with 'val'
  SparseRowMatrix(int rows, int cols, int nz, const T& val);
 
  //! reserve memory for given rows, and number of non zeros, 
  //! set all entries to value of val
  void reserve(int rows, int cols, int nz, const T& val);

  //! resize keeping old values if possible, assuming rows == cols  
  void resize ( int newSize );

  //! resize keeping old values if possible   
  void resize ( int newRow, int newCol );

  //! free memory for values_ and col_
  ~SparseRowMatrix();
    
/*******************************/
/*  Access and info functions  */
/*******************************/
  //! length of array used to store matrix 
  int numberOfValues () const { return dim_[0]*nz_; }
  //! matrix value taken from real array 
  T&      val(int i) 
  { 
    assert( i >= 0 );
    assert( i <= numberOfValues () );
    return values_[i]; 
  }

  //! return value and clear matrix entry 
  T popValue (int i) { 
    T v = val(i); 
    values_[i] = 0;
    col_[i] = -1;
    return v;
  }

  //! return real column number for (row,localCol) 
  int realCol (int row, int fakeCol) const 
  {
    int pos = row*nz_ + fakeCol;
    return col_[pos];
  }

  //! return pair< value, column >, used by BlockMatrix 
  std::pair < T , int > realValue(int row, int fakeCol) 
  {
    assert( fakeCol < nz_ );
    int pos = row*nz_ + fakeCol;
    return realValue(pos);
  }
  
  //! return pair< value, column >, used by BlockMatrix 
  std::pair < const T , int > realValue(int row, int fakeCol) const 
  {
    assert( fakeCol < nz_ );
    int pos = row*nz_ + fakeCol;
    return realValue(pos);
  }
  
  //! return pair< value, column >, used by BlockMatrix 
  std::pair < T , int > realValue(int index) 
  {
    return std::pair< T, int > (values_[index], col_[index]); 
  }

  //! return pair< value, column >, used by BlockMatrix 
  std::pair < const T , int > realValue(int index) const 
  {
    return std::pair< const T, int > (values_[index], col_[index]); 
  }

  //! returns local col index for given global (row,col) 
  int colIndex(int row, int col);
  
  //! returns true if entry (row,col) exists in matrix 
  bool find (int row, int col) const;

  //! return reference to value on given entry 
  const T&  val(int i) const { return values_[i]; }

  //! return number of rows = 0, cols = 1
  int dim(int i) const {return dim_[i];}
  //! return number of rows = 0, cols = 1
  int size(int i) const {return dim_[i];}
  
  //! return number of rows  
  int rows() const {return dim_[0];}

  //! return number of columns 
  int cols() const {return dim_[1];}

  //! return max number of non zeros 
  int NumNonZeros() const {return nz_;}
  
  //! return number of non zeros in row 
  int NumNonZeros(int i) const 
  { 
    assert( nonZeros_ );
    return nonZeros_[i]; 
  }

  //! return value of entry (i,j)
  T  operator() (int i, int j) const;        

  //! set entry to value 
  void set(int row, int col, T val);
  
  //! set all entries in row to zero 
  void clearRow (int row);
  //! set all matrix entries to zero 
  void clear();
  
  //! add value to row,col entry 
  void add(int row, int col, T val);
  //! muliply with scalar value 
  void multScalar(int row, int col, T val);
  
  //! make unitRow(row) and unitCol(col)
  void kroneckerKill(int row, int col);

  //! same as apply A * x = ret 
  template <class VECtype> 
  void mult(const VECtype *x, VECtype * ret) const;

  //! same as apply A * x = ret, used by OEM-Solvers 
  template <class VECtype> 
  void multOEM(const VECtype *x, VECtype * ret) const;

  //! same as apply A^T * x = ret, used by OEM-Solvers 
  template <class VECtype> 
  void multOEM_t(const VECtype *x, VECtype * ret) const;

  //! A(f) = ret, same as mult 
  template <class DiscFType, class DiscFuncType>
  void apply(const DiscFType &f, DiscFuncType &ret) const;
  
  //! A^T(f) = ret
  template <class DiscFuncType>
  void apply_t(const DiscFuncType &f, DiscFuncType &ret) const;
  
  //! A(f) = ret 
  template <class DiscFuncType> 
  void operator () (const DiscFuncType &f, DiscFuncType &ret) const 
  {
    apply(f,ret); 
  };
  
  //! return diagonal of (this * A * B)
  template <class DiscFuncType>
  void getDiag(const ThisType&A, const ThisType &B, DiscFuncType &rhs) const;
  
  //! return diagonal of (this * A)
  template <class DiscFuncType>
  void getDiag(const ThisType &A, DiscFuncType &rhs) const;
  
  //! return diagonal entries of this matrix 
  template <class DiscFuncType>
  void getDiag(DiscFuncType &rhs) const;
  
  //! add diagonal to given DiscreteFunction
  template <class DiscFuncType>
  void addDiag(DiscFuncType &rhs) const;
  
  //! print matrix 
  void print (std::ostream& s) const;

  //! print values 
  void printReal (std::ostream& s) const;
     
  //! make row a row with 1 on diagonal and all other entries 0 
  void unitRow(int row);

  //! make column a column with 1 on diagonal and all other entries 0 
  void unitCol(int col);

  //! check symetry 
  void checkSym ();

  // res = this * B 
  void multiply(const ThisType & B, ThisType & res) const;
  //! add other matrix to this matrix 
  void add(const ThisType & B ); 

  //! resort to have ascending column numbering 
  void resort();
  //! resort row to have ascending column numbering 
  void resortRow(const int row);

private:
  //! delete memory 
  void removeObj();
};

} // end namespace Sparselib

#include "spmatrix.cc"

#endif  
