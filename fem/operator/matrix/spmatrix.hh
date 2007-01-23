#ifndef DUNE_SPMATRIX_HH
#define DUNE_SPMATRIX_HH

#include <vector>

#include <dune/fem/discretefunction/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/space/common/communicationmanager.hh>

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

  // omega for ssor preconditioning
  const double omega_;
  
  //! Copy Constructor prohibited 
  SparseRowMatrix(const SparseRowMatrix<T> &S);
public:
  //! makes Matrix of zero length, omega is 1.1 by default 
  SparseRowMatrix(double omega = 1.1); 

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
  T&  val(int i) 
  { 
    assert( i >= 0 );
    assert( i < numberOfValues () );
    return values_[i]; 
  }

  //! return reference to value on given entry 
  const T&  val(int i) const { 
    assert( i >= 0 );
    assert( i < numberOfValues () );
    return values_[i]; 
  }

  //! return value and clear matrix entry 
  T popValue (int i) 
  { 
    T v = val(i); 
    values_[i] = 0;
    col_[i] = -1;
    return v;
  }

  //! return real column number for (row,localCol) 
  int realCol (int row, int fakeCol) const 
  {
    assert( row < dim_[0] );
    int pos = row*nz_ + fakeCol;
    return col_[pos];
  }

  //! return pair< value, column >, used by BlockMatrix 
  std::pair < T , int > realValue(int row, int fakeCol) 
  {
    assert( row < dim_[0] );
    assert( fakeCol < nz_ );
    int pos = row*nz_ + fakeCol;
    return realValue(pos);
  }
  
  //! return pair< value, column >, used by BlockMatrix 
  std::pair < const T , int > realValue(int row, int fakeCol) const 
  {
    assert( row < dim_[0] );
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
  T multOEMRow (const VECtype *x , const int row ) const;

  //! same as apply A * x = ret, used by OEM-Solvers 
  template <class VECtype> 
  void multOEM(const VECtype *x, VECtype * ret) const;

  //! calculates ret += A * x 
  template <class VECtype> 
  void multOEMAdd(const VECtype *x, VECtype * ret) const;

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

  //! SSOR preconditioning 
  void ssorPrecondition (const T*, T*) const;

  //! returns true if preconditioing is called before matrix multiply 
  bool rightPrecondition() const { return true; }
  
  //! apply preconditioning, calls ssorPreconditioning at the moment 
  void precondition (const T*u , T*x) const {  ssorPrecondition(u,x); }

private:
  //! delete memory 
  void removeObj();
};

template <class RowSpaceImp, class ColumnSpaceImp> 
class SparseRowMatrixObject
{
public:
  typedef RowSpaceImp RowSpaceType;
  typedef ColumnSpaceImp ColumnSpaceType;

private:
  typedef typename RowSpaceType::GridType::template Codim<0>::Entity EntityType;
  typedef SparseRowMatrixObject<RowSpaceType,ColumnSpaceType> ThisType;
public:  
  typedef SparseRowMatrix<double> MatrixType;
  typedef MatrixType PreconditionMatrixType;
  //typedef DFAdapt<RowSpaceType> PreconditionMatrixType;
  
  //! LocalMatrix 
  template <class MatrixImp> 
  class LocalMatrix
  {
    typedef MatrixImp MatrixType;
    MatrixType & matrix_; 
    
    const int rowIndex_;
    const int colIndex_;

    std::vector<int> row_;
    std::vector<int> col_;
    
  public:  
    //! constructor taking entity and spaces for using mapToGlobal
      //, class RowSpaceType, class ColSpaceType> 
    LocalMatrix(MatrixType & m,
                const EntityType& rowEntity,
                const RowSpaceType & rowSpace,
                const EntityType& colEntity,
                const ColumnSpaceType & colSpace)
      : matrix_(m)
      , rowIndex_(rowSpace.indexSet().index(rowEntity))
      , colIndex_(colSpace.indexSet().index(colEntity))
    {
      row_.resize(rowSpace.getBaseFunctionSet(rowEntity).numBaseFunctions());
      col_.resize(colSpace.getBaseFunctionSet(colEntity).numBaseFunctions());

      {
        const size_t rows = row_.size();
        for(size_t i=0; i<rows; ++i) 
          row_[i] = rowSpace.mapToGlobal( rowEntity, i );
      }
      {
        const size_t cols = col_.size();
        for(size_t i=0; i<cols; ++i) 
          col_[i] = colSpace.mapToGlobal( colEntity, i );
      }
    }

  private: 
    //! copy not allowed 
    LocalMatrix(const LocalMatrix &);

  public:
    //! return number of rows 
    int rows () const { return row_.size(); }
    //! return number of cols 
    int cols () const { return col_.size(); }

    //! add value to matrix entry
    void add(int localRow, int localCol , const double value)
    {
      assert( localRow >= 0 );
      assert( localCol >= 0 );

      assert( localRow < (int) row_.size() );
      assert( localCol < (int) col_.size() );
      matrix_.add(row_[localRow],col_[localCol],value);
    }

    //! get matrix entry 
    double get(int localRow, int localCol) const 
    {
      assert( localRow >= 0 );
      assert( localCol >= 0 );

      assert( localRow < (int) row_.size() );
      assert( localCol < (int) col_.size() );
      return matrix_(row_[localRow],col_[localCol]);
    }
    
    //! set matrix enrty to value 
    void set(int localRow, int localCol, const double value)
    {
      assert( localRow >= 0 );
      assert( localCol >= 0 );

      assert( localRow < (int) row_.size() );
      assert( localCol < (int) col_.size() );
      matrix_.set(row_[localRow],col_[localCol],value);
    }

    //! clear all entries belonging to local matrix 
    void clear ()
    {
      const int row = rows();
      for(int i=0; i<row; ++i)
      {
        matrix_.clearRow( row_[i] );
      }
    }

    //! resort all global rows of matrix to have ascending numbering 
    void resort ()
    {
      const int row = rows();
      for(int i=0; i<row; ++i)
      {
        matrix_.resortRow( row_[i] );
      }
    }
  };

public:
  typedef LocalMatrix<MatrixType> LocalMatrixType;

// commented out, as AdaptiveDiscreteFunction is not known as type
  typedef AdaptiveDiscreteFunction<RowSpaceType> DestinationType;

  typedef CommunicationManager<RowSpaceType> CommunicationManagerType; 

  const RowSpaceType & rowSpace_; 
  const ColumnSpaceType & colSpace_;
  
  int rowMaxNumbers_;
  int sequence_;

  MatrixType matrix_; 
  bool preconditioning_;
  PreconditionMatrixType * pcMatrix_;

  mutable CommunicationManagerType communicate_;

  //! setup matrix handler 
  SparseRowMatrixObject(const RowSpaceType & rowSpace, 
                        const ColumnSpaceType & colSpace,
                        const std::string& paramfile ) 
    : rowSpace_(rowSpace)
    , colSpace_(colSpace) 
    , rowMaxNumbers_(-1)
    , sequence_(-1)
    , matrix_()
    , preconditioning_(false)
    , pcMatrix_(0)
    , communicate_(rowSpace_)
  {
    if( paramfile != "" )
    {
      int precon = 0;
      readParameter(paramfile,"Preconditioning",precon);
      preconditioning_ = (precon == 1) ? true : false;
    }
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

  PreconditionMatrixType& pcMatrix () { 
    //assert( pcMatrix_ );
    //return *pcMatrix_; 
    return matrix_;
  }

  //! reserve memory corresponnding to size of spaces 
  void reserve(bool verbose = false ) 
  {
    if(sequence_ != rowSpace_.sequence())
    {
      // if empty grid do nothing (can appear in parallel runs)
      if( (rowSpace_.begin() != rowSpace_.end()) && 
          (colSpace_.begin() != colSpace_.end()) )
      {
        
        rowMaxNumbers_ = rowSpace_.getBaseFunctionSet(*(rowSpace_.begin())).numBaseFunctions();

        if(verbose) 
        {
          std::cout << "Reserve Matrix with (" << rowSpace_.size() << "," << colSpace_.size()<< ")\n";
          std::cout << "Number of base functions = (" << rowMaxNumbers_ << ")\n";
        }

        assert( rowMaxNumbers_ > 0 );

        // factor for non-conforming grid is 4 in 3d and 2 in 2d  
        //const int factor = (Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));
        const int factor = 1; //(Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));

        // upper estimate for number of neighbors 
        enum { dim = RowSpaceType :: GridType :: dimension };
        rowMaxNumbers_ *= (factor * dim * 2) + 1; // e.g. 7 for dim = 3

        matrix_.reserve(rowSpace_.size(),colSpace_.size(),rowMaxNumbers_,0.0);

        /*
        if(hasPcMatrix())
        {
          pcMatrix_ = new PreconditionMatrixType("pcMatrix",rowSpace_);
        }
        */
      }
      sequence_ = rowSpace_.sequence();
    }
  }

  //! mult method of matrix object used by oem solver
  void multOEM(const double * arg, double * dest) const 
  {
    communicate( arg );
    matrix_.multOEM(arg,dest);
  }

  //! communicate data 
  void communicate(const double * arg) const
  {
    if( rowSpace_.grid().comm().size() <= 1 ) return ;

    DestinationType tmp("SparseRowMatrixObject::communicate_tmp",rowSpace_,arg);
    communicate_.exchange( tmp );
  }


  //! resort row numbering in matrix to have ascending numbering 
  void resort() 
  {
    matrix_.resort();
  }

  void createPreconditionMatrix()
  { 
    /*
    if(hasPcMatrix())
    {
      PreconditionMatrixType & diag = pcMatrix(); 
      diag.clear();
      
      matrix_.addDiag( diag );
  
      double * diagPtr = diag.leakPointer();
      const int singleSize = rowSpace_.size();
      for(register int i=0; i<singleSize; ++i) 
      {
        double val = diagPtr[i];
        // when using parallel Version , we could have zero on diagonal
        // for ghost elements 
        //assert( (spc_.grid().comm().size() > 1) ? 1 : (std::abs( val ) > 0.0
        if( std::abs( val ) > 0.0 )
        {
          val = 1.0/val;
          diagPtr[i] = val;
        }
        else
          diagPtr[i] = 1.0;
      }
    }
    */
  }

};

} // end namespace Dune 

#include "spmatrix.cc"
#endif  
