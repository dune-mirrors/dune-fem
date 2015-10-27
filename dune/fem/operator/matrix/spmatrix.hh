#ifndef DUNE_FEM_SPMATRIX_HH
#define DUNE_FEM_SPMATRIX_HH

//- system includes
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>

//- local includes
#include <dune/fem/function/adaptivefunction/adaptivefunction.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/io/file/asciiparser.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/storage/objectstack.hh>

namespace Dune
{

  namespace Fem
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
    //! Compressed row sparse matrix, where only the nonzeros of a row are keeped
    //! (except if you "set" a single element explicitly with the value 0)
    //!
    //*****************************************************************

    template <class T>
    class SparseRowMatrix
    {
      enum { defaultCol = -1 };
      enum { firstCol = defaultCol + 1 };

    public:
      //! matrix field type
      typedef T Ttype;
      typedef SparseRowMatrix<T> ThisType;
      //! type of the base matrix
      //! for consistency with ISTLMatrixObject
      typedef ThisType MatrixBaseType;

      SparseRowMatrix(const SparseRowMatrix<T> & ) = delete;

      //! make matrix of zero length
      explicit SparseRowMatrix(double omega = 1.1);

      //! make matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      //! and intialize all values with 'val'
      SparseRowMatrix(int rows, int cols, int nz, const T& val = 0, double omega = 1.1 );

      //! free memory
      ~SparseRowMatrix();

      //! reserve memory for given rows, and number of non zeros
      void reserve(int rows, int cols, int nz,const T& dummy);

      //! return number of rows
      int rows() const
      {
        return dim_[0];
      }

      //! return number of columns
      int cols() const
      {
        return dim_[1];
      }

      //! set entry to value (also setting 0 will result in an entry)
      void set(int row, int col, T val);

      //! add value to row,col entry
      void add(int row, int col, T val);

      //! A(f) = ret, same as mult
      template <class DiscFType, class DiscFuncType>
      void apply(const DiscFType &f, DiscFuncType &ret) const;

      //! return value of entry (row,col)
      T operator()( const int row, const int col ) const;
      T operator()( const unsigned int row, const unsigned int col ) const
      {
        return (*this)( int( row ), int( col ) );
      }
      T operator()( const long unsigned int row, const long unsigned int col ) const
      {
        return this->operator()((unsigned int)(row), (unsigned int)(col) );
      }

      //! set all entries in row to zero
      void clearRow (int row);

      //! set all matrix entries to zero
      void clear();

      //! return max number of non zeros
      //! used in SparseRowMatrixObject::reserve
      int numNonZeros() const {return nz_;}

      //! return number of non zeros in row
      //! used in ColCompMatrix::setMatrix
      int numNonZeros(int i) const
      {
        assert( nonZeros_ );
        return nonZeros_[i];
      }

      //! return pair< value, column >
      //! used in ColCompMatrix::setMatrix
      std::pair < const T , int > realValue(int index) const
      {
        return std::pair< const T, int > (values_[index], col_[index]);
      }

    private:
      //! resize keeping old values if possible
      void resize( int newRow, int newCol, int newNz );

      //! returns local col index for given global (row,col)
      int colIndex(int row, int col);

      //! check whether number of stored nonzeros corresponds to the counters in nonZeros[]
      //! and check, whether all columns within this range have reasonable values
      //! return true if consistent and false if not consistent
      //! an assert(checkConsistency()) can be called at entry and exit of
      //! non-const sparsematrix operations for ensuring maintaining of
      //! consistence. This can be made conditional by the member variable
      //! checkNonConstMethods
      bool checkConsistency() const;

      //! return real column number for (row,localCol)
      int realCol(int row, int fakeCol) const
      {
        assert(fakeCol<dim_[1]);
        assert( row < dim_[0] );
        int pos = row*nz_ + fakeCol;
        return col_[pos];
      }

      //! make row a row with 1 on diagonal and all other entries 0
      void unitRow(int row);

      //! make column a column with 1 on diagonal and all other entries 0
      void unitCol(int col);

      //! check symetry
      void checkSym();

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
      void ssorPrecondition(const T* , T* ) const;

      //! returns true if preconditioing is called before matrix multiply
      bool rightPrecondition() const
      {
        return true;
      }

      //! apply preconditioning, calls ssorPreconditioning at the moment
      void precondition(const T *u , T *x) const
      {
        ssorPrecondition(u,x);
      }

    private:
      //! delete memory
      void removeObj();

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
      double omega_;
    };

    template< class DomainSpace, class RangeSpace,
              class Matrix = SparseRowMatrix< typename DomainSpace :: RangeFieldType > >
    class SparseRowMatrixObject
    {
    public:
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      /*******************************************************************
      *   Rows belong to the DomainSpace and Columns to the RangeSpace   *
      *******************************************************************/
      typedef typename DomainSpaceType :: EntityType  DomainEntityType ;
      typedef typename RangeSpaceType  :: EntityType  RangeEntityType ;
      typedef typename DomainSpaceType :: EntityType  ColumnEntityType ;
      typedef typename RangeSpaceType  :: EntityType  RowEntityType ;

      typedef typename DomainSpaceType :: BlockMapperType DomainBlockMapperType ;
      typedef NonBlockMapper< DomainBlockMapperType, DomainSpaceType :: localBlockSize > DomainMapperType;
      typedef typename RangeSpaceType :: BlockMapperType RangeBlockMapperType ;
      typedef NonBlockMapper< RangeBlockMapperType, RangeSpaceType :: localBlockSize > RangeMapperType;
      typedef Matrix MatrixType ;
      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, MatrixType > ThisType;

    protected:
      typedef typename DomainSpaceType :: GridType GridType;

      template< class MatrixObject >
      struct LocalMatrixTraits;

      template< class MatrixObject >
      class LocalMatrix;

    public:
      typedef MatrixType PreconditionMatrixType;

      //! type of local matrix
      typedef LocalMatrix<ThisType> ObjectType;
      typedef ThisType LocalMatrixFactoryType;
      typedef Fem :: ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
      //! type of local matrix
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

      typedef ColumnObject< ThisType > LocalColumnObjectType;

    protected:
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;

      DomainMapperType domainMapper_ ;
      RangeMapperType  rangeMapper_ ;

      int sequence_;

      mutable MatrixType matrix_;
      bool preconditioning_;

      mutable LocalMatrixStackType localMatrixStack_;

    public:
      //! setup matrix handler
      inline SparseRowMatrixObject( const DomainSpaceType &domainSpace,
                                    const RangeSpaceType &rangeSpace,
                                    const std::string &paramfile = "" )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        domainMapper_( domainSpace_.blockMapper() ),
        rangeMapper_( rangeSpace_.blockMapper() ),
        sequence_( -1 ),
        matrix_(),
        preconditioning_( false ),
        localMatrixStack_( *this )
      {
        int precon = 0;
        if( paramfile != "" )
        {
          readParameter( paramfile, "Preconditioning", precon );
        }
        else
        {
          precon = Parameter :: getValue("Preconditioning", precon );
        }
        preconditioning_ = (precon > 0) ? true : false;
      }

      //! return domain space (i.e. space that builds the rows)
      const DomainSpaceType& domainSpace() const
      {
        return domainSpace_;
      }

      //! return range space (i.e. space that builds the columns)
      const RangeSpaceType& rangeSpace() const
      {
        return rangeSpace_;
      }

      //! return reference to storage object
      MatrixType &matrix() const
      {
        return matrix_;
      }

      //! interface method from LocalMatrixFactory
      inline ObjectType *newObject() const
      {
        return new ObjectType( *this, domainSpace_, rangeSpace_, domainMapper_, rangeMapper_ );
      }

      //! return local matrix
      inline LocalMatrixType localMatrix( const DomainEntityType &domainEntity,
                                          const RangeEntityType &rangeEntity ) const
      {
        return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
      }

      LocalColumnObjectType localColumn( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType( *this, domainEntity );
      }

      //! resize all matrices and clear them
      inline void clear()
      {
        matrix_.clear();
      }

      //! reserve memory for assemble based on the provided stencil
      template <class Stencil>
      inline void reserve(const Stencil &stencil, bool verbose = false )
      {
        if( sequence_ != domainSpace_.sequence() )
        {
#ifndef DUNE_FEM_DONT_CHECKITERATORS_OF_SPACE
          // if empty grid do nothing (can appear in parallel runs)
          if( (domainSpace_.begin() != domainSpace_.end())
              && (rangeSpace_.begin() != rangeSpace_.end()) )
#endif
          {

            if( verbose )
            {
              const int rowMaxNumbers = rangeMapper_.maxNumDofs();
              const int colMaxNumbers = domainMapper_.maxNumDofs();

              std::cout << "Reserve Matrix with (" << rangeSpace_.size() << "," << domainSpace_.size()<< ")\n";
              std::cout << "Max number of base functions = (" << rowMaxNumbers << "," << colMaxNumbers << ")\n";
            }

            // upper estimate for number of non-zeros
            const static size_t domainLocalBlockSize = DomainSpaceType::localBlockSize;
            const int nonZeros = std::max( (int)(stencil.maxNonZerosEstimate()*domainLocalBlockSize),
                                           matrix_.numNonZeros() );
            matrix_.reserve( rangeSpace_.size(), domainSpace_.size(), nonZeros, 0.0 );
          }
          sequence_ = domainSpace_.sequence();
        }
      }

      //! apply matrix to discrete function
      template< class DomainFunction, class RangeFunction >
      void apply( const DomainFunction &arg, RangeFunction &dest ) const
      {
        // do matrix vector multiplication
        matrix_.apply( arg, dest );

        // communicate data
        dest.communicate();
      }

      ///////////////////////////////////////////////////
      // Method used for the diagonal preconditioning of the dune fem
      // inbuild solvers
      ///////////////////////////////////////////////////
      //! extract diagonal entries from matrix into discrete function
      template < class DiscreteFunctionType >
      void extractDiagonal( DiscreteFunctionType& diag ) const
      {
        // this only works for matrices with same number of rows,cols
        assert( matrix_.rows() == matrix_.cols() );
        typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType ;
        const DofIteratorType dofEnd = diag.dend();
        unsigned int row = 0;
        for( DofIteratorType dofIt = diag.dbegin();
             dofIt != dofEnd; ++ dofIt, ++row )
        {
          assert( row < ( unsigned int )matrix_.rows() );
          (*dofIt) = matrix_( row, row );
        }
      }

      //! resort row numbering in matrix to have ascending numbering
      void resort()
      {
        matrix_.resort();
      }

      //! mult method of matrix object used by oem solver
      double ddotOEM( const double *v, const double *w ) const
      {
        typedef AdaptiveDiscreteFunction< DomainSpace > DomainFunctionType;
        DomainFunctionType V( "ddot V", domainSpace_, v );
        DomainFunctionType W( "ddot W", domainSpace_, w );
        return V.scalarProductDofs( W );
      }

      //! mult method of matrix object used by oem solver
      void multOEM( const double *arg, double *dest ) const
      {
        typedef AdaptiveDiscreteFunction< DomainSpace > DomainFunctionType;
        typedef AdaptiveDiscreteFunction< RangeSpace > RangeFunctionType;

        DomainFunctionType farg( "multOEM arg", domainSpace_, arg );
        RangeFunctionType fdest( "multOEM dest", rangeSpace_, dest );
        apply( farg, fdest );
      }
    };



    template< class DomainSpace, class RangeSpace, class Matrix >
    template< class MatrixObject >
    struct SparseRowMatrixObject< DomainSpace, RangeSpace, Matrix >::LocalMatrixTraits
    {
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, Matrix > SparseRowMatrixObjectType;

      typedef typename SparseRowMatrixObjectType :: template LocalMatrix< MatrixObject > LocalMatrixType;

      typedef typename RangeSpaceType :: RangeFieldType RangeFieldType;
      typedef RangeFieldType LittleBlockType;

      typedef typename SparseRowMatrixObjectType::DomainMapperType  DomainMapperType;
      typedef typename SparseRowMatrixObjectType::RangeMapperType   RangeMapperType;
    };



    //! LocalMatrix
    template< class DomainSpace, class RangeSpace, class Matrix >
    template< class MatrixObject >
    class SparseRowMatrixObject< DomainSpace, RangeSpace, Matrix > :: LocalMatrix
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

      //! type of nonblocked domain mapper
      typedef typename Traits :: DomainMapperType DomainMapperType;
      //! type of nonblocked domain mapper
      typedef typename Traits :: RangeMapperType RangeMapperType;

    protected:
      MatrixType &matrix_;
      const DomainMapperType& domainMapper_;
      const RangeMapperType&  rangeMapper_;

      typedef std :: vector< typename RangeMapperType :: SizeType > RowIndicesType ;
      //! global index in the DomainSpace
      RowIndicesType rowIndices_;

      typedef std :: vector< typename DomainMapperType :: SizeType > ColumnIndicesType ;
      //! global index in the RangeSpace
      ColumnIndicesType columnIndices_;

      using BaseType :: domainSpace_;
      using BaseType :: rangeSpace_;

    public:
      //! constructor
      inline LocalMatrix( const MatrixObjectType &matrixObject,
                          const DomainSpaceType &domainSpace,
                          const RangeSpaceType &rangeSpace,
                          const DomainMapperType& domainMapper,
                          const RangeMapperType& rangeMapper )
      : BaseType( domainSpace, rangeSpace),
        matrix_( matrixObject.matrix() ),
        domainMapper_( domainMapper ),
        rangeMapper_( rangeMapper )
      {
      }

      LocalMatrix( const LocalMatrix & ) = delete;

      void init( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        // initialize base functions sets
        BaseType::init( domainEntity, rangeEntity );

        // rows are determined by the range space
        rowIndices_.resize( rangeMapper_.numDofs( rangeEntity ) );
        rangeMapper_.mapEach( rangeEntity, AssignFunctor< RowIndicesType >( rowIndices_ ) );

        // columns are determind by the domain space
        columnIndices_.resize( domainMapper_.numDofs( domainEntity ) );
        domainMapper_.mapEach( domainEntity, AssignFunctor< ColumnIndicesType >( columnIndices_ ) );
      }

      //! return number of rows
      int rows() const
      {
        return rowIndices_.size();
      }

      //! return number of columns
      int columns() const
      {
        return columnIndices_.size();
      }

      //! add value to matrix entry
      void add( int localRow, int localCol, const DofType value )
      {
        assert( value == value );
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        matrix_.add( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! get matrix entry
      DofType get( int localRow, int localCol ) const
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        return matrix_( rowIndices_[ localRow ], columnIndices_[ localCol ] );
      }

      //! set matrix entry to value
      void set( int localRow, int localCol, const DofType value )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        matrix_.set( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! set matrix row to zero except diagonla entry
      void unitRow( const int localRow )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.unitRow( rowIndices_[ localRow ] );
      }

      //! set matrix row to zero
      void clearRow( const int localRow )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.clearRow( rowIndices_[localRow]);
      }

      //! set matrix column to zero
      void clearCol( const int localCol )
      {
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.clearCol( columnIndices_[localCol] );
      }

      //! clear all entries belonging to local matrix
      void clear()
      {
        const int row = rows();
        for( int i = 0; i < row; ++i )
          matrix_.clearRow( rowIndices_[ i ] );
      }

      //! scale local matrix with a certain value
      void scale( const DofType& value )
      {
        const int row = rows();
        for( int i = 0; i < row; ++i )
          matrix_.scaleRow( rowIndices_[ i ] , value );
      }

      //! resort all global rows of matrix to have ascending numbering
      void resort()
      {
        const int row = rows();
        for( int i = 0; i < row; ++i )
          matrix_.resortRow( rowIndices_[ i ] );
      }

    };

//////////////////////////////////////////////////////////////////////////////////////////////////
/// Out of class defined methods for the SparseMatrix class
//////////////////////////////////////////////////////////////////////////////////////////////////

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
      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    template <class T>
    SparseRowMatrix<T>::SparseRowMatrix(int rows, int cols, int nz,
                                        const T& dummy, double omega)
    : omega_(omega)
    {
      // standard settings as above
      values_ = 0;
      col_ = 0;
      dim_[0] = 0;
      dim_[1] = 0;
      memSize_ = 0;
      nz_ = 0;
      nonZeros_ = 0;

      // resize and get storage
      reserve(rows,cols,nz,dummy);

      // set all values to default value
      clear();

      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::removeObj()
    {
      if(checkNonConstMethods)
        assert(checkConsistency());
      if(values_)
        delete [] values_;
      if(col_)
        delete [] col_;
      if(nonZeros_)
        delete [] nonZeros_;
      values_ = 0;
      col_ = 0;
      nonZeros_ = 0;
      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    template <class T>
    SparseRowMatrix<T>::~SparseRowMatrix()
    {
      if(checkNonConstMethods)
        assert(checkConsistency());
      removeObj();
      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    /***********************************/
    /*  Construct from storage vectors */
    /***********************************/
    template <class T>
    void SparseRowMatrix<T>::reserve(int rows, int cols, int nz,const T& dummy)
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

      memSize_ = rows * nz;
      nz_ = nz;
      // add first col for offset
      nz_ += firstCol ;

      assert( dim_[0] > 0 );
      assert( dim_[1] > 0 );

      // make resize
      newValues_.resize( nz_ );

      // only reserve for indices
      newIndices_.reserve( nz_ );

      // set all values to default value
      clear();

      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    // resize matrix
    template <class T>
    void SparseRowMatrix<T>::resize(int newRow, int newCol, int newNz )
    {
      if(newRow != this->rows() || newNz > nz_ )
      {
        if( newNz < 0 )
          newNz = nz_;

        int newMemSize = newRow * newNz ;

        int memHalf = (int) memSize_/2;
        if((newMemSize > memSize_) || (newMemSize < memHalf))
        {
          T tmp = 0;
          T * oldValues = values_;       values_ = 0;
          int * oldCol  = col_;          col_ = 0;
          int * oldNonZeros = nonZeros_; nonZeros_ = 0;
          const int oldNz = nz_;
          const int copySize = std::min( dim_[0] , newRow );
          const int oldSize = dim_[0];

          // reserve new memory
          reserve(newRow,newCol,newNz,tmp);

          if( (oldSize > 0) && (oldNz > 0 ))
          {
            std::memset( col_ , -1 , newRow * newNz * sizeof(int));
            for( int row = 0; row < copySize; ++ row )
            {
              const int newLoc = row * newNz ;
              const int oldLoc = row * oldNz ;
              std::memcpy( values_ + newLoc , oldValues + oldLoc , oldNz * sizeof(T) );
              std::memcpy( col_ + newLoc    , oldCol + oldLoc   , oldNz * sizeof(int) );
            }
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

      assert( this->rows()  == newRow );
      assert( this->cols()  == newCol );
    }

    template< class T >
    inline T SparseRowMatrix<T>::operator()( const int row, const int col ) const
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
      if(checkNonConstMethods)
        assert(checkConsistency());
      assert( row >= 0 );
      assert( row < dim_[0] );

      int i = 0;
      while ( i < nz_ && col_[row*nz_+i] < col && col_[row*nz_+i] != defaultCol )
        ++i;
      if(col_[row*nz_+i] == col)
        return i;  // column already in matrix
      else if( col_[row*nz_+i] == defaultCol )
      { // add this column at end of this row
        ++nonZeros_[row];
        return i;
      }
      else
      {
        ++nonZeros_[row];
        // must shift this row to add col at the position i
        int j = nz_-1; // last column
        if (col_[row*nz_+j] != defaultCol)
        { // new space available - so resize
          resize( rows(), cols(), (2 * nz_) );
          j++;
        }
        for(;j>i;--j)
        {
          col_[row*nz_+j] = col_[row*nz_+j-1];
          values_[row*nz_+j] = values_[row*nz_+j-1];
        }
        col_[row*nz_+i] = col;
        values_[row*nz_+i] = 0;
        return i;
      }
    }

    template <class T>
    void SparseRowMatrix<T>::clear()
    {
      T init = 0;
      for(int i=0; i<dim_[0]*nz_; ++i)
      {
        values_ [i] = init;
        col_[i] = defaultCol;
      }

      for(int i=0; i<dim_[0]; ++i)
      {
        nonZeros_[i] = 0;
      }

      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::clearRow(int row)
    {
      if(checkNonConstMethods)
        assert(checkConsistency());
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

      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::set(int row, int col, T val)
    {
      if(checkNonConstMethods)
        assert(checkConsistency());
      assert((col>=0) && (col <= dim_[1]));
      assert((row>=0) && (row <= dim_[0]));

      int whichCol = colIndex(row,col);
      assert( whichCol != defaultCol );

      values_[row*nz_ + whichCol] = val;
      if(whichCol >= nonZeros_[row])
        nonZeros_[row]++;
      col_[row*nz_ + whichCol] = col;

      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    template <class T>
    void SparseRowMatrix<T>::add(int row, int col, T val)
    {
      if (checkNonConstMethods) assert(checkConsistency());
      int whichCol = colIndex(row,col);
      assert( whichCol != defaultCol );
      values_[row*nz_ + whichCol] += val;
      col_[row*nz_ + whichCol] = col;
      if(checkNonConstMethods)
        assert(checkConsistency());
    }

    /***************************************/
    /*  Matrix-MV_Vector multiplication    */
    /***************************************/
    template <class T> template <class ArgDFType, class DestDFType>
    void SparseRowMatrix<T>::apply(const ArgDFType &f, DestDFType &ret) const
    {
      typedef typename DestDFType::DofIteratorType DofIteratorType;

      typedef typename ArgDFType :: ConstDofBlockPtrType ConstDofBlockPtrType;
      enum { blockSize = ArgDFType :: DiscreteFunctionSpaceType :: localBlockSize };

      //! we assume that the dimension of the functionspace of f is the same as
      //! the size of the matrix
      DofIteratorType ret_it = ret.dbegin();

      for(int row=0; row<dim_[0]; ++row)
      {
        (*ret_it) = 0.0;

        //! DofIteratorType schould be the same
        for(int col=firstCol; col<nz_; ++col)
        {
          const int thisCol = row*nz_ + col;
          const int realCol = col_[thisCol];

          if( realCol == defaultCol )
            continue;

          const int blockNr = realCol / blockSize ;
          const int dofNr = realCol % blockSize ;
          ConstDofBlockPtrType fBlock = f.block( blockNr );
          (*ret_it) += values_[thisCol] * (*fBlock)[ dofNr ];
        }

        ++ret_it;
      }
    }

    template <class T>
    bool SparseRowMatrix<T>::checkConsistency() const
    {
      bool consistent = true;

      // only perform check, if there is any data:
      if(nonZeros_ || values_ || col_)
      {
        for(int row=0; row<dim_[0]; row++)
        {
          if(nonZeros_[row]<0 || nonZeros_[row]> dim_[1])
          {
            std::cout << "error in consistency of row " << row
              << ": NonZeros_[row] = "<< nonZeros_[row]
              << " is not within reasonable range "
              << " of dim = (" << dim_[0]<<","<< dim_[1]<< ")"
              << std::endl;
            consistent = false;
            return consistent;
          }

          for (int fakeCol =0; fakeCol < nonZeros_[row]; fakeCol++)
            if ((realCol(row,fakeCol)<0) || (realCol(row,fakeCol)>=dim_[1]))
            {
              std::cout << "error in consistency of row " << row
                << ": NonZeros_[row] = "<< nonZeros_[row]
                << ", fakeCol = " << fakeCol << ", realCol(fakeCol) = "
                << realCol(row,fakeCol) << std::endl;
              consistent = false;
              return consistent;
            }
        }
        return consistent;
      }

      assert(consistent);

      return consistent;
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

          if(realCol < row)
            dot += localValues[col] * x[realCol];
          else if(realCol == row)
          {
            diag = localValues[col];
            assert( std::abs(diag) > 0.0 );
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

          if(realCol > row)
            dot += localValues[col] * x[realCol];
          else if(realCol == row)
          {
            diag = localValues[col];
            assert( std::abs(diag) > 0.0 );
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

          if(realCol > row)
            dot += localValues[col] * x[realCol];
          else if(realCol == row)
          {
            diag = localValues[col];
            assert( std::abs(diag) > 0.0 );
          }
          ++thisCol;
        }
        x[row] -= omega * dot / diag;
      }
    }
  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPMATRIX_HH
