#ifndef DUNE_ONTHEFLYMATRIX_HH
#define DUNE_ONTHEFLYMATRIX_HH

//- system includes 
#include <stack>

//- local includes 
#include "blockmatrix.hh"

namespace Dune
{

//*****************************************************************
//
//  --OnTheFlyMatrix
//    
//*****************************************************************
template <class T>
class OnTheFlyMatrix 
{
public: 
  typedef T Ttype;  //! remember the value type
  typedef OnTheFlyMatrix<T> ThisType;

private:
  typedef DenseMatrix<T> DenseMatrixType;
  typedef SparseRowMatrix< DenseMatrixType * > OnTheFlyMatrixType;

  typedef std::stack< DenseMatrixType * > MatrixStackType;

  OnTheFlyMatrixType matrix_;
  int localRows_;
  int localCols_;

  bool rowWise_; 

  mutable MatrixStackType freeStack_;

  int nonZeros_;

  OnTheFlyMatrix(const OnTheFlyMatrix<T> &S); //! Copy Constructor
public:
  //! makes Matrix of zero length
  OnTheFlyMatrix() {}

  //! make matrix with 'rows' rows and 'cols' columns,
  //! maximum 'nz' non zero values in each row 
  //! and intialize all values with 'val'
  
  void reserve(int rows, int cols, int lrows, int lcols, 
               int nonZeros, bool rowWise = true) 
  { 
  }
  
  void clear () 
  {
  }
    
/*******************************/
/*  Access and info functions  */
/*******************************/
  int size(int i) const { return 0; }

  int littleRows() const { return localRows_;}
  int littleCols() const { return localCols_;}

  void set(int row, int col, DenseMatrixType & val)
  {
  }
  
  void add(int row, int col, DenseMatrixType & val)
  {
  }

  void resize( int nsize) 
  {
  }

  void resize( int rows, int cols ) 
  {
  }

  void resort () 
  {
  }

  void resortRow ( int row ) 
  {
  }
  
  void multiply(const ThisType & A , const ThisType & B )
  {
  }

  void add(const ThisType & A)
  {
  }

  // result = this * vec 
  void multOEM(const T * vec, T * result) const
  {
    abort();
  }
  
  // result = this * vec 
  void multOEMAdd(const T * vec, T * result) const
  {
    abort();
  }
  
  //! this precondition is from right side 
  bool rightPrecondition () const { return true; }
  
  // result = this * vec 
  void precondition(const T * u, T * x) const
  {
  }
  
  // result = this * vec 
  void setDiag(const T * diag)
  {
  }

  void getDiag(T * diag) 
  {
  }
  
  void print(std::ostream & s) const 
  {
  }

};

template <class RowSpaceImp, class ColumnSpaceImp> 
class OnTheFlyMatrixObject;

template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
struct OnTheFlyMatrixTraits
{
  typedef RowSpaceImp RowSpaceType;
  typedef ColSpaceImp ColumnSpaceType;
  typedef OnTheFlyMatrixTraits<RowSpaceType,ColumnSpaceType> ThisType;

  template <class OperatorTraits>
  struct MatrixObject
  {
    typedef OnTheFlyMatrixObject<RowSpaceType,ColumnSpaceType> MatrixObjectType;
  };
};



//! matrix object holding a blockamtrix
template <class RowSpaceImp, class ColumnSpaceImp> 
class OnTheFlyMatrixObject
{
public:  
  typedef RowSpaceImp RowSpaceType;
  typedef ColumnSpaceImp ColumnSpaceType ;

  //! number of rows of blocks 
  enum { littleRows = RowSpaceType :: localBlockSize };
  //! number of columns of blocks 
  enum { littleCols = ColumnSpaceType :: localBlockSize };

  typedef typename RowSpaceType::GridType::template Codim<0>::Entity EntityType;

  typedef OnTheFlyMatrixObject<RowSpaceType,ColumnSpaceType> ThisType;

  typedef OnTheFlyMatrix<double> MatrixType;
  typedef MatrixType PreconditionMatrixType;

  template <class MatrixObjectImp> 
  class LocalMatrix; 

  struct LocalMatrixTraits
  {
    typedef RowSpaceImp DomainSpaceType ;
    typedef ColumnSpaceImp RangeSpaceType;
    typedef typename RowSpaceImp :: RangeFieldType RangeFieldType;
    typedef LocalMatrix<ThisType>  LocalMatrixType;
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
    typedef typename RowSpaceType :: RangeFieldType DofType;

    //! type of little blocks 
    typedef DenseMatrix<DofType> LittleMatrixType;

    //! type of row mapper 
    typedef typename MatrixObjectType :: RowMapperType RowMapperType;
    //! type of col mapper 
    typedef typename MatrixObjectType :: ColMapperType ColMapperType;

  private:  
    const MatrixObjectType & matrixObj_; 
    LittleMatrixType matrix_; 

  public:  
    LocalMatrix(const MatrixObjectType & mObj,
                const RowSpaceType & rowSpace,
                const ColumnSpaceType & colSpace)
      : BaseType(rowSpace,colSpace) 
      , matrixObj_(mObj)
      , matrix_() 
    {
    }

    //! initialize local matrix to entities 
    void init(const EntityType& rowEntity,
              const EntityType& colEntity)
    {
      // initialize base functions sets 
      BaseType :: init ( rowEntity , colEntity );

      const int numRows = this->rows(); 
      const int numCols = this->columns();

      matrix_.resize( numRows, numCols );
    }

  private: 
    //! prohibit copy cosntructor 
    LocalMatrix(const LocalMatrix &);

    // check whether given (row,col) pair is valid
    void check(int localRow, int localCol) const
    {
      assert( localRow >= 0 && localRow < matrix_.rows());
      assert( localCol >= 0 && localCol < matrix_.cols());
    }
    
  public:
    //! add value to matrix 
    void add(int localRow, int localCol , const DofType value)
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      matrix_[localRow][localCol] += value;
    }

    //! return matrix entry 
    const DofType get(int localRow, int localCol ) const 
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      return matrix_[localRow][localCol];
    }

    //! set matrix enrty to value 
    void set(int localRow, int localCol, const DofType value)
    {
#ifndef NDEBUG
      check(localRow,localCol);
#endif
      matrix_[localRow][localCol] = value;
    }

    //! set matrix enrty to value 
    void unitRow(const int localRow) 
    {
      const int col = this->columns();
      for(int localCol=0; localCol<col; ++localCol)
      {
        matrix_[localRow][localCol] = 0;
      }
      matrix_[localRow][localRow] = 1;
    }

    //! clear all entries belonging to local matrix 
    void clear ()
    {
      matrix_.clear();
    }

    //! resort does nothing because already sorted in good way 
    void resort ()
    {
    }
  };

public:
  //! type of local matrix 
  typedef LocalMatrix<ThisType> ObjectType;
  typedef ThisType LocalMatrixFactoryType;
  typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
  //! type of local matrix 
  typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

  typedef typename RowSpaceType :: BlockMapperType RowMapperType; 
  typedef typename ColumnSpaceType :: BlockMapperType ColMapperType; 

  typedef AdaptiveDiscreteFunction<RowSpaceType> DestinationType;

  typedef CommunicationManager<RowSpaceType> CommunicationManagerType;

  // row space 
  const RowSpaceType & rowSpace_; 
  // column space 
  const ColumnSpaceType & colSpace_;

  mutable MatrixType matrix_;

  // communication manager 
  mutable CommunicationManagerType communicate_;

  // local matrix stack 
  mutable LocalMatrixStackType localMatrixStack_;

  //! constructor 
  //! \param rowSpace space defining row structure 
  //! \param colSpace space defining column structure 
  //! \param paramfile parameter file to read variables 
  //!         - Preconditioning: {0 == no, 1 == yes} used is SSOR 
  OnTheFlyMatrixObject(const RowSpaceType & rowSpace, 
                       const ColumnSpaceType & colSpace,
                       const std::string& paramfile) 
    : rowSpace_(rowSpace)
    , colSpace_(colSpace) 
    , matrix_()
    , communicate_(rowSpace_)
    , localMatrixStack_(*this)
  {
    /*
    if( paramfile != "" )
    {
      int precon = 0;
      readParameter(paramfile,"Preconditioning",precon);
      preconditioning_ = (precon > 0) ? true : false;
    } 
    */
    assert( rowSpace_.indexSet().size(0) == 
            colSpace_.indexSet().size(0) ); 
  }

  //! return reference to stability matrix 
  MatrixType & matrix() const 
  { 
    return matrix_; 
  }

  //! set all matrix entries to zero  
  void clear() 
  {
  }

  //! return true if precoditioning matrix is provided 
  bool hasPreconditionMatrix () const { return false; }

  //! return reference to preconditioner (here also systemMatrix)
  const PreconditionMatrixType& preconditionMatrix () const { return matrix(); }

  //! reserve memory corresponnding to size of spaces 
  void reserve(bool verbose = false ) 
  {
  }

  //! mult method of matrix object used by oem solver
  void multOEM(const double * arg, double * dest) const 
  {
    abort();
  }

  //! communicate data 
  void communicate(const double * arg) const
  {
    /*
    if( rowSpace_.grid().comm().size() <= 1 ) return ;

    DestinationType tmp("BlockMatrixObject::communicate_tmp",rowSpace_,arg);
    communicate_.exchange( tmp );
    */
  }

  //! resort row numbering in matrix to have ascending numbering 
  void resort() 
  {
  }

  //! print matrix 
  void print(std::ostream & s) const 
  {
  }

 //! interface method from LocalMatrixFactory 
  ObjectType* newObject() const
  {
    return new ObjectType(*this,
                          rowSpace_,
                          colSpace_);
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
