#ifndef DUNE_SPMATRIX_HH
#define DUNE_SPMATRIX_HH

//- system includes 
#include <vector>

//- local includes 
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/operator/common/localmatrix.hh> 
#include <dune/fem/operator/common/localmatrixwrapper.hh> 

namespace Dune
{
  
  //! anonymous namespace, such that variable is only known within 
  //! SparseRowMatrix and other classes located here
  namespace
  {    
//! If you have problems with this sparsematrix class, this might be due to
//! inconsistencies produced in some methods.  
//! In this case, you should turn on the consistencycheck of all non-const 
//! methods by setting the following variable to 1 / 0 and by this locate the 
//! buggy member method. Default is 0 = Check Off
    const int checkNonConstMethods = 0;
  }
  
//*****************************************************************
//
//  --SparseRowMatrix
//  
//! Compressed row sparse matrix, where only the nonzeros of a row are
//! keeped 
//! (except if you "set" a single element explicitly 
//! with the value 0, which is not forbidden and an element entry is 
//! created)
//!
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
  T* values_ ;      //! data values (nz_ * dim_[0] elements)
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
  //! and intialize all values with 'val' and set omega_ to omega
  SparseRowMatrix(int rows, int cols, int nz, const T& val = 0, 
                  double omega = 1.1);
 
  //! reserve memory for given rows, and number of non zeros, 
  //! set all entries to value dummy.... What is the use of this value?
  //! only initializing with 0 makes sense, so this argument is renamed
  //! 'dummy', by the way, nothing is happening with this argument, so 
  //! might be completely removed.
  void reserve(int rows, int cols, int nz, const T& dummy);

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
    if (checkNonConstMethods) assert(checkConsistency());
    T v = val(i); 
    values_[i] = 0;
    col_[i] = -1;
    if (checkNonConstMethods) assert(checkConsistency());
    return v;
  }

  //! return real column number for (row,localCol) 
  int realCol (int row, int fakeCol) const 
  {
    assert(fakeCol<dim_[1]);
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
  int numNonZeros() const {return nz_;}
  
  //! return number of non zeros in row 
  int numNonZeros(int i) const 
  { 
    assert( nonZeros_ );
    return nonZeros_[i]; 
  }

  //! return value of entry (i,j)
  T  operator() (int i, int j) const;        

  //! set entry to value
  //! note, that every entry is performed into the matrix!
  //! also setting of value 0 will result in an entry. So these
  //! calls should be ommited on a higher level 
  void set(int row, int col, T val);
  
  //! set all entries in row to zero 
  void clearRow (int row);

  //! set all matrix entries to zero, no other value makes sense for 
  //! sparse matrix 
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
     
  //! print columns
  void printColumns(std::ostream& s) const;
  
  //! print row-wise stored number of nonzeros 
  //! No counting is performed, but the member variable nonZeros_[]
  //! is reported. So here inconsistencies can occur to the true 
  //! nonzero entries in the matrix.
  void printNonZeros(std::ostream& s) const;
  
  //! check consistency, i.e. whether number of stored nonzeros
  //! corresponds to the counters in nonZeros[] and check, whether all
  //! columns within this range have reasonable values
  //! true == consistent
  //! false == non consistent
  //! an assert(checkConsistency()) can be called at entry and exit of
  //! non-const sparsematrix operations for ensuring maintaining of
  //! consistence. This can be made conditional by the member variable 
  //! checkNonConstMethods
  bool checkConsistency() const;

  //! make row a row with 1 on diagonal and all other entries 0 
  void unitRow(int row);

  //! make column a column with 1 on diagonal and all other entries 0 
  void unitCol(int col);

  //! check symetry 
  void checkSym ();

  // res = this * B 
  void multiply(const ThisType & B, ThisType & res) const;

  //! multiply this matrix with scalar 
  void scale(const T& factor);
  
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

  template <class DomainSpace, class RangeSpace>
  class SparseRowMatrixObject;

  template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
  struct SparseRowMatrixTraits
  {
    typedef RowSpaceImp RowSpaceType;
    typedef ColSpaceImp ColumnSpaceType;
    typedef SparseRowMatrixTraits<RowSpaceType,ColumnSpaceType> ThisType;

    template <class OperatorTraits>
    struct MatrixObject
    {
      typedef SparseRowMatrixObject<RowSpaceType,ColumnSpaceType> MatrixObjectType;
    };
  };

  template< class DomainSpace, class RangeSpace >
  class SparseRowMatrixObject
  {
  public:
    typedef DomainSpace DomainSpaceType;
    typedef RangeSpace RangeSpaceType;

  private:  
    typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType > ThisType;

  protected:
    typedef typename DomainSpaceType :: GridType GridType;
    typedef typename GridType :: template Codim< 0 > :: Entity EntityType;

    template< class MatrixObject >
    struct LocalMatrixTraits;

    template< class MatrixObject >
    class LocalMatrix;
    
  public:  
    typedef SparseRowMatrix< double > MatrixType;
    typedef MatrixType PreconditionMatrixType;

  public:
    //! type of local matrix 
    typedef LocalMatrix<ThisType> ObjectType;
    typedef ThisType LocalMatrixFactoryType;
    typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
    //! type of local matrix 
    typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

    typedef CommunicationManager< RangeSpaceType > CommunicationManagerType;

  protected:
    const DomainSpaceType &domainSpace_;
    const RangeSpaceType &rangeSpace_;
    
    int rowMaxNumbers_;
    int sequence_;

    mutable MatrixType matrix_;
    bool preconditioning_;
    PreconditionMatrixType * pcMatrix_;

    mutable CommunicationManagerType communicate_;

    mutable LocalMatrixStackType localMatrixStack_;

  public:
    //! setup matrix handler 
    inline SparseRowMatrixObject( const DomainSpaceType &domainSpace,
                                  const RangeSpaceType &rangeSpace,
                                  const std :: string &paramfile = "" )
    : domainSpace_( domainSpace ),
      rangeSpace_( rangeSpace ),
      rowMaxNumbers_( -1 ),
      sequence_( -1 ),
      matrix_(),
      preconditioning_( false ),
      pcMatrix_( 0 ),
      communicate_( rangeSpace_ , InteriorBorder_All_Interface ),
      localMatrixStack_( *this )
    {
      if( paramfile != "" )
      {
        int precon = 0;
        readParameter( paramfile, "Preconditioning", precon );
        preconditioning_ = (precon > 0 ? true : false);
      }
    }

    //! return reference to stability matrix
    inline MatrixType &matrix () const
    {
      return matrix_;
    }

    //! interface method from LocalMatrixFactory
    inline ObjectType *newObject () const
    {
      return new ObjectType( *this, domainSpace_, rangeSpace_ );
    }

    //! return local matrix 
    inline LocalMatrixType localMatrix( const EntityType &rowEntity,
                                        const EntityType &colEntity ) const
    {
      return LocalMatrixType( localMatrixStack_, rowEntity, colEntity );
    }

    //! resize all matrices and clear them 
    inline void clear ()
    {
      matrix_.clear();
    }

    //! return true if precoditioning matrix is provided 
    bool hasPreconditionMatrix () const
    {
      return preconditioning_;
    }

    //! return reference to preconditioner 
    const PreconditionMatrixType &preconditionMatrix () const 
    { 
      return matrix_;
    }

    //! reserve memory corresponnding to size of spaces
    inline void reserve(bool verbose = false )
    {
      if( sequence_ != domainSpace_.sequence() )
      {
        // if empty grid do nothing (can appear in parallel runs)
        if( (domainSpace_.begin() != domainSpace_.end())
            && (rangeSpace_.begin() != rangeSpace_.end()) )
        {
          
          rowMaxNumbers_ = domainSpace_.baseFunctionSet(*(domainSpace_.begin())).numBaseFunctions();
          int colMaxNumbers = rangeSpace_.baseFunctionSet(*(domainSpace_.begin())).numBaseFunctions();

          rowMaxNumbers_ = std::max(rowMaxNumbers_, colMaxNumbers);

          if( verbose )
          {
            std::cout << "Reserve Matrix with (" << domainSpace_.size() << "," << rangeSpace_.size()<< ")\n";
            std::cout << "Number of base functions = (" << rowMaxNumbers_ << ")\n";
          }

          assert( rowMaxNumbers_ > 0 );

          // factor for non-conforming grid is 4 in 3d and 2 in 2d  
          //const int factor = (Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));
          const int factor = 1; //(Capabilities::isLeafwiseConforming<GridType>::v) ? 1 : (2 * (dim-1));

          // upper estimate for number of neighbors 
          enum { dim = GridType :: dimension };
          rowMaxNumbers_ *= (factor * dim * 2) + 1; // e.g. 7 for dim = 3

          matrix_.reserve( domainSpace_.size(), rangeSpace_.size(), rowMaxNumbers_, 0.0 );
        }
        sequence_ = domainSpace_.sequence();
      }
    }

    //! apply matrix to discrete function
    template< class DomainFunction, class RangeFunction >
    void apply ( const DomainFunction &arg, RangeFunction &dest ) const
    {
      matrix_.multOEM( arg.leakPointer(), dest.leakPointer() );
      communicate_.exchange( dest, (DFCommunicationOperation :: Add *) 0 );
    }

    //! mult method of matrix object used by oem solver
    void multOEM( const double *arg, double *dest ) const
    {
      typedef AdaptiveDiscreteFunction< DomainSpaceType > DomainFunctionType;
      typedef AdaptiveDiscreteFunction< RangeSpaceType > RangeFunctionType;

      DomainFunctionType farg( "multOEM arg", domainSpace_, arg );
      RangeFunctionType fdest( "multOEM dest", rangeSpace_, dest );
      apply( farg, fdest );
    }

    //! resort row numbering in matrix to have ascending numbering 
    void resort() 
    {
      matrix_.resort();
    }

    void createPreconditionMatrix()
    { 
    }
  };



  template< class DomainSpace, class RangeSpace >
  template< class MatrixObject >
  struct SparseRowMatrixObject< DomainSpace, RangeSpace > :: LocalMatrixTraits
  {
    typedef DomainSpace DomainSpaceType;
    typedef RangeSpace RangeSpaceType;
    
    typedef typename SparseRowMatrixObject< DomainSpaceType, RangeSpaceType >
      :: template LocalMatrix< MatrixObject > LocalMatrixType;

    typedef typename RangeSpaceType :: RangeFieldType RangeFieldType;
    typedef RangeFieldType LittleBlockType;
  };



  //! LocalMatrix 
  template< class DomainSpace, class RangeSpace >
  template< class MatrixObject >
  class SparseRowMatrixObject< DomainSpace, RangeSpace > :: LocalMatrix
  : public LocalMatrixDefault< LocalMatrixTraits< MatrixObject > >
  {
  public:
    //! type of matrix object 
    typedef MatrixObject MatrixObjectType;

    //! type of the traits
    typedef LocalMatrixTraits< MatrixObjectType > Traits;

  private:
    typedef LocalMatrixDefault< Traits > BaseType;

  public:
    //! type of matrix 
    typedef typename MatrixObjectType :: MatrixType MatrixType;

    //! type of entries of little blocks 
    typedef typename Traits :: RangeFieldType RangeFieldType;

    //! type of the DoFs
    typedef RangeFieldType DofType;

    //! type of little blocks 
    typedef typename Traits :: LittleBlockType LittleBlockType;

  protected:
    MatrixType &matrix_; 
    
    //! global row numbers 
    std :: vector< int > row_;
    //! global col numbers  
    std :: vector< int > col_;

    using BaseType :: domainSpace_;
    using BaseType :: rangeSpace_;
    
  public:  
    //! constructor taking entity and spaces for using mapToGlobal
    //! class RowSpaceType, class ColSpaceType> 
    inline LocalMatrix( const MatrixObjectType &matrixObject,
                        const DomainSpaceType &domainSpace,
                        const RangeSpaceType &rangeSpace )
    : BaseType( domainSpace, rangeSpace),
      matrix_( matrixObject.matrix() )
    {
    }
    
  private: 
    // prohibit copying 
    LocalMatrix( const LocalMatrix & );

  public:
    void init( const EntityType &rowEntity, const EntityType &colEntity )
    {
      // initialize base functions sets 
      BaseType :: init ( rowEntity , colEntity );
        
      row_.resize( domainSpace_.baseFunctionSet( rowEntity ).numBaseFunctions() );
      col_.resize( rangeSpace_.baseFunctionSet( colEntity ).numBaseFunctions() );

      const size_t rows = row_.size();
      for( size_t i = 0; i < rows; ++i )
        row_[ i ] = domainSpace_.mapToGlobal( rowEntity, i );
      
      const size_t cols = col_.size();
      for( size_t i = 0; i < cols; ++i )
        col_[ i ] = rangeSpace_.mapToGlobal( colEntity, i );
    }

    //! return number of rows 
    int rows () const
    {
      return row_.size();
    }

    //! return number of columns 
    int columns () const
    {
      return col_.size();
    }

    //! add value to matrix entry
    void add( int localRow, int localCol, const DofType value )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );

      matrix_.add( row_[ localRow ], col_[ localCol ], value );
    }

    //! get matrix entry 
    DofType get( int localRow, int localCol ) const
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );

      return matrix_( row_[ localRow ], col_[ localCol ] );
    }
    
    //! set matrix entry to value 
    void set( int localRow, int localCol, const DofType value )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );

      matrix_.set( row_[ localRow ], col_[ localCol ], value );
    }

    //! set matrix row to zero except diagonla entry 
    void unitRow( const int localRow )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      matrix_.unitRow( row_[ localRow ] );
    }

    //! clear all entries belonging to local matrix 
    void clear ()
    {
      const int row = rows();
      for( int i = 0; i < row; ++i )
        matrix_.clearRow( row_[ i ] );
    }

    //! resort all global rows of matrix to have ascending numbering 
    void resort ()
    {
      const int row = rows();
      for( int i = 0; i < row; ++i )
        matrix_.resortRow( row_[ i ] );
    }
  };

} // end namespace Dune 

#include "spmatrix.cc"

#endif  
