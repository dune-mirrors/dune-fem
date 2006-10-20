#ifndef DUNE_SPMATRIX_HH
#define DUNE_SPMATRIX_HH

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

private:
  T* values_ ;       //! data values (nz_ elements)
  int* col_;        //! row_ptr (dim_[0]+1 elements)
  int dim_[2];      //! dim_[0] x dim_[1] Matrix
  int nz_;          //! number of nonzeros per row

  int memSize_;
  
public:
  SparseRowMatrix(); //! makes Matrix of zero length
  SparseRowMatrix(const SparseRowMatrix<T> &S); //! Copy Constructor

  //! make matrix with 'rows' rows and 'cols' columns,
  //! maximum 'nz' non zero values in each row 
  //! and intialize all values with 'val'
  SparseRowMatrix(int rows, int cols, int nz, const T& val);
  
  void makeSpMat(int rows, int cols, int nz, const T& val);

  //! free memory for values_ and col_
  ~SparseRowMatrix();
    
/*******************************/
/*  Access and info functions  */
/*******************************/
  // length of array used to store matrix 
  int numberOfValues () const { return dim_[0]*nz_; }
  // matrix value taken from real array 
  T&      val(int i) 
  { 
    assert( i >= 0 );
    assert( i <= numberOfValues () );
    return values_[i]; 
  }

  // return value and clear matrix entry 
  T popValue (int i) { 
    T v = val(i); 
    values_[i] = 0;
    col_[i] = -1;
    return v;
  }

  int realCol (int row, int fakeCol) const 
  {
    int pos = row*nz_ + fakeCol;
    return col_[pos];
  }

  std::pair < T , int > realValue(int row, int fakeCol) 
  {
    assert( fakeCol < nz_ );
    int pos = row*nz_ + fakeCol;
    return realValue(pos);
  }
  
  std::pair < const T , int > realValue(int row, int fakeCol) const 
  {
    assert( fakeCol < nz_ );
    int pos = row*nz_ + fakeCol;
    return realValue(pos);
  }
  
  std::pair < T , int > realValue(int index) 
  {
    return std::pair< T, int > (values_[index], col_[index]); 
  }

  std::pair < const T , int > realValue(int index) const 
  {
    return std::pair< const T, int > (values_[index], col_[index]); 
  }

  int colIndex(int row, int col);
  
  bool find (int row, int col) const;

  const T&  val(int i) const { return values_[i]; }

  int dim(int i) const {return dim_[i];};
  int size(int i) const {return dim_[i];};
  int NumNonZeros() const {return nz_;};

  T  operator() (int i, int j) const;        

  void set(int row, int col, T val);
  void setLastRow(int col, T val);
  
  void clear();
  
  void add(int row, int col, T val);
  void multScalar(int row, int col, T val);
   
  
  void kroneckerKill(int row, int col);

  // same as apply A * x = ret 
  template <class VECtype> 
  void mult(const VECtype *x, VECtype * ret) const;

  // same as apply A * x = ret 
  template <class VECtype> 
  void multOEM(const VECtype *x, VECtype * ret) const;

  // same as apply A * x = ret 
  template <class VECtype> 
  void multOEM_t(const VECtype *x, VECtype * ret) const;

  // A f = ret, same as mult 
  template <class DiscFType, class DiscFuncType>
  void apply(const DiscFType &f, DiscFuncType &ret) const;
  
  template <class DiscFuncType>
  void apply_t(const DiscFuncType &f, DiscFuncType &ret) const;
  
  template <class DiscFuncType>
  void prepareGlobal(const DiscFuncType &f, DiscFuncType &ret) const {}
  
  void finalizeGlobal() const {}
  
  template <class DiscFuncType> 
  void operator () (const DiscFuncType &f, DiscFuncType &ret) const 
  {
    apply(f,ret); 
  };
  
  template <class DiscFuncType>
  void diagCond (DiscFuncType &rhs);
  
  template <class DiscFuncType>
  void getDiag(const ThisType&, const ThisType &, DiscFuncType &rhs) const;
  
  template <class DiscFuncType>
  void getDiag(const ThisType &, DiscFuncType &rhs) const;
  
  template <class DiscFuncType>
  void getDiag(DiscFuncType &rhs) const;
  
  //! add diagonal to given DiscreteFunction
  template <class DiscFuncType>
  void addDiag(DiscFuncType &rhs) const;
  
  void print (std::ostream& s) const;
  void printReal (std::ostream& s) const;
     
  void unitRow(int row);
  void unitCol(int col);

  void checkSym ();

  void resize ( int newSize );

  void resize ( int newRow, int newCol );

  // res = this * B 
  void multiply(const ThisType & B, ThisType & res) const;
  void add(const ThisType & B ); 

private:
   
};

} // end namespace Sparselib

#include "spmatrix.cc"

#endif  
