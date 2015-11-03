#ifndef DUNE_FEM_SPMATRIX_HH
#define DUNE_FEM_SPMATRIX_HH

// C++ includes
#include <vector>
#include <algorithm>
#include <iostream>
#include <array>
#include <string>

// DUNE-FEM includes
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

    //! SparseRowMatrix
    template <class T>
    class SparseRowMatrix
    {
      static constexpr int defaultCol = -1;
      static constexpr int firstCol = defaultCol + 1;

    public:
      //! matrix field type
      typedef T Ttype;
      typedef SparseRowMatrix<T> ThisType;
      //! type of the base matrix
      //! for consistency with ISTLMatrixObject
      typedef ThisType MatrixBaseType;

      SparseRowMatrix(const ThisType& ) = delete;

      //! construct matrix of zero size
      explicit SparseRowMatrix() :
        values_(0), col_(0), nonZeros_(0), dim_({{0,0}}), nz_(0)
      {}

      //! construct matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      SparseRowMatrix(int rows, int cols, int nz) :
      values_(0), col_(0), nonZeros_(0), dim_({{0,0}}), nz_(0)
      {
        reserve(rows,cols,nz);
      }

      //! reserve memory for given rows, columns and number of non zeros
      void reserve(int rows, int cols, int nz, T&& )
      {
        if( (rows != dim_[0]) || (cols != dim_[1]) || (nz != nz_))
          resize(rows,cols,nz);
        clear();
      }

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
      void set(int row, int col, T val)
      {
        assert((col>=0) && (col <= dim_[1]));
        assert((row>=0) && (row <= dim_[0]));

        auto whichCol = colIndex(row,col);
        assert( whichCol != defaultCol );

        values_[row*nz_ + whichCol] = val;
        if(whichCol >= nonZeros_[row])
          nonZeros_[row]++;
        col_[row*nz_ + whichCol] = col;
      }

      //! add value to row,col entry
      void add(int row, int col, T val)
      {
        assert((col>=0) && (col <= dim_[1]));
        assert((row>=0) && (row <= dim_[0]));

        auto whichCol = colIndex(row,col);
        assert( whichCol != defaultCol );

        values_[row*nz_ + whichCol] += val;
        col_[row*nz_ + whichCol] = col;
      }

      //! ret = A*f
      template<class ArgDFType, class DestDFType>
      void apply(const ArgDFType& f, DestDFType& ret) const
      {
        constexpr auto blockSize = ArgDFType :: DiscreteFunctionSpaceType :: localBlockSize;
        auto ret_it = ret.dbegin();

        for(auto row=0; row<dim_[0]; ++row)
        {
          (*ret_it) = 0.0;
          for(auto col=firstCol; col<nz_; ++col)
          {
            const auto thisCol = row*nz_ + col;
            const auto realCol = col_[thisCol];

            if( realCol == defaultCol )
              continue;

            const auto blockNr = realCol / blockSize ;
            const auto dofNr = realCol % blockSize ;
            auto fBlock = f.block( blockNr );
            (*ret_it) += values_[thisCol] * (*fBlock)[ dofNr ];
          }

          ++ret_it;
        }
      }

      //! return value of entry (row,col)
      T operator()(int row, int col) const
      {
        assert((col>=0) && (col <= dim_[1]));
        assert((row>=0) && (row <= dim_[0]));

        const auto nonZ = nonZeros_[row];
        auto thisCol = row*nz_;
        for(auto i=firstCol; i<nonZ; ++i)
        {
          if(col_[thisCol] == col)
            return values_[thisCol];
          ++thisCol;
        }
        return 0;
      }

      //! set all matrix entries to zero
      void clear()
      {
        for(auto& entry : values_)
          entry = 0;
        for(auto& entry : col_)
          entry = defaultCol;
        for(auto& entry : nonZeros_)
          entry = 0;
      }

      //! set all entries in row to zero
      void clearRow(int row)
      {
        assert((row>=0) && (row <= dim_[0]));

        nonZeros_[row] = firstCol;
        auto col = row * nz_;
        for(auto i=0; i<nz_; ++i)
        {
          values_ [col] = 0;
          col_[col] = defaultCol;
          ++col;
        }
      }

      //! return max number of non zeros
      //! used in SparseRowMatrixObject::reserve
      int numNonZeros() const
      {
        return nz_;
      }

      //! return number of non zeros in row
      //! used in ColCompMatrix::setMatrix
      int numNonZeros(int i) const
      {
        return nonZeros_[i];
      }

      //! return pair (value,column)
      //! used in ColCompMatrix::setMatrix
      std::pair<const T, int> realValue(int index) const
      {
        return std::pair<const T, int>(values_[index], col_[index]);
      }

    private:
      //! resize matrix
      void resize(int rows, int cols, int nz)
      {
        constexpr auto colVal(defaultCol);
        values_.resize( rows*nz , 0 );
        col_.resize( rows*nz , colVal );
        nonZeros_.resize( rows , 0 );
        dim_[0] = rows;
        dim_[1] = cols;
        nz_ = nz+firstCol;
      }

      //! returns local col index for given global (row,col)
      int colIndex(int row, int col)
      {
        assert((col>=0) && (col <= dim_[1]));
        assert((row>=0) && (row <= dim_[0]));

        int i = 0;
        while( i < nz_ && col_[row*nz_+i] < col && col_[row*nz_+i] != defaultCol )
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
          auto j = nz_-1; // last column
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

      //! return real column number for (row,localCol)
      int realCol(int row, int fakeCol) const
      {
        assert( fakeCol < dim_[1] );
        assert( row < dim_[0] );
        auto pos = row*nz_ + fakeCol;
        return col_[pos];
      }

      std::vector<T> values_;
      std::vector<int> col_;
      std::vector<int> nonZeros_;
      std::array<int,2> dim_;
      int nz_;
    };



    //! SparseRowMatrixObject
    template< class DomainSpace, class RangeSpace,
              class Matrix = SparseRowMatrix< typename DomainSpace :: RangeFieldType > >
    class SparseRowMatrixObject
    {
    public:
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;
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
          readParameter( paramfile, "Preconditioning", precon );
        else
          precon = Parameter :: getValue("Preconditioning", precon );
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
              const auto rowMaxNumbers = rangeMapper_.maxNumDofs();
              const auto colMaxNumbers = domainMapper_.maxNumDofs();

              std::cout << "Reserve Matrix with (" << rangeSpace_.size() << "," << domainSpace_.size()<< ")" << std::endl;
              std::cout << "Max number of base functions = (" << rowMaxNumbers << "," << colMaxNumbers << ")" << std::endl;
            }

            // upper estimate for number of non-zeros
            constexpr auto domainLocalBlockSize = DomainSpaceType::localBlockSize;
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

      //! extract diagonal entries from matrix into discrete function
      template < class DiscreteFunctionType >
      void extractDiagonal( DiscreteFunctionType& diag ) const
      {
        // this only works for matrices with same number of rows,cols
        assert( matrix_.rows() == matrix_.cols() );
        const auto dofEnd = diag.dend();
        int row = 0;
        for( auto dofIt = diag.dbegin(); dofIt != dofEnd; ++dofIt, ++row )
        {
          assert( row < matrix_.rows() );
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



    //! LocalMatrixTraits
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
      {}

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
      void add(int localRow, int localCol, DofType value)
      {
        assert( value == value );
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        matrix_.add( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! get matrix entry
      DofType get(int localRow, int localCol) const
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        return matrix_( rowIndices_[ localRow ], columnIndices_[ localCol ] );
      }

      //! set matrix entry to value
      void set(int localRow, int localCol, DofType value)
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );

        matrix_.set( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! set matrix row to zero except diagonla entry
      void unitRow(int localRow)
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.unitRow( rowIndices_[ localRow ] );
      }

      //! set matrix row to zero
      void clearRow( int localRow )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.clearRow( rowIndices_[localRow]);
      }

      //! set matrix column to zero
      void clearCol( int localCol )
      {
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.clearCol( columnIndices_[localCol] );
      }

      //! clear all entries belonging to local matrix
      void clear()
      {
        const auto row = rows();
        for( auto i = 0; i < row; ++i )
          matrix_.clearRow( rowIndices_[ i ] );
      }

      //! scale local matrix with a certain value
      void scale( const DofType& value )
      {
        const auto row = rows();
        for( auto i = 0; i < row; ++i )
          matrix_.scaleRow( rowIndices_[ i ] , value );
      }

      //! resort all global rows of matrix to have ascending numbering
      void resort()
      {
        const auto row = rows();
        for( auto i = 0; i < row; ++i )
          matrix_.resortRow( rowIndices_[ i ] );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPMATRIX_HH
