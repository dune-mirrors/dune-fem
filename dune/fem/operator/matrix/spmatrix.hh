#ifndef DUNE_FEM_SPMATRIX_HH
#define DUNE_FEM_SPMATRIX_HH

// system includes
#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

// local includes
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

#include <dune/fem/operator/common/stencil.hh>

#include <dune/fem/operator/matrix/functor.hh>

namespace Dune
{

  namespace Fem
  {

    struct MatrixParameter
      : public Dune::Fem::LocalParameter< MatrixParameter, MatrixParameter >
    {

      MatrixParameter( const std::string keyPrefix = "" )
      {}

      MatrixParameter( const ParameterReader &parameter, const std::string keyPrefix = "istl." )
      {}

      virtual double overflowFraction () const
      {
        return 1.0;
      }

      virtual int numIterations () const
      {
        return 0;
      }

      virtual double relaxation () const
      {
        return 1.0;
      }

      virtual int method () const
      {
        return 0;
      }

      virtual bool fastILUStorage () const { return false; }

      virtual std::string preconditionName() const
      {
        return "None";
      }

      virtual bool viennaCL () const { return false; }
      virtual bool blockedMode () const { return false; }
    };


    struct SparseRowMatrixParameter
      : public MatrixParameter
    {
      typedef MatrixParameter BaseType;

      SparseRowMatrixParameter( const std::string keyPrefix = "spmatrix." )
        : BaseType( keyPrefix )
      {}
    };




    //! SparseRowMatrix
    template <class T>
    class SparseRowMatrix
    {
    public:
      //! matrix field type
      typedef T field_type;
      //! matrix index type
      typedef std::size_t size_type;
      typedef SparseRowMatrix<field_type> ThisType;
      //! type of the base matrix
      //! for consistency with ISTLMatrixObject
      typedef ThisType MatrixBaseType;

      static constexpr size_type defaultCol = std::numeric_limits<size_type>::max();
      static constexpr size_type firstCol = 0;

      SparseRowMatrix(const ThisType& ) = delete;

      //! construct matrix of zero size
      explicit SparseRowMatrix() :
        values_(0), col_(0), nonZeros_(0), dim_({{0,0}}), nz_(0)
      {}

      //! construct matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      SparseRowMatrix(size_type rows, size_type cols, size_type nz) :
        values_(0), col_(0), nonZeros_(0), dim_({{0,0}}), nz_(0)
      {
        reserve(rows,cols,nz);
      }

      //! reserve memory for given rows, columns and number of non zeros
      void reserve(size_type rows, size_type cols, size_type nz)
      {
        if( (rows != dim_[0]) || (cols != dim_[1]) || (nz != nz_))
          resize(rows,cols,nz);
        clear();
      }

      //! return number of rows
      size_type rows() const
      {
        return dim_[0];
      }

      //! return number of columns
      size_type cols() const
      {
        return dim_[1];
      }

      //! set entry to value (also setting 0 will result in an entry)
      void set(size_type row, size_type col, field_type val)
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
      void add(size_type row, size_type col, field_type val)
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
        constexpr auto blockSize = ArgDFType::DiscreteFunctionSpaceType::localBlockSize;
        auto ret_it = ret.dbegin();

        for(auto row=decltype(dim_[0]){0}; row<dim_[0]; ++row)
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
            (*ret_it) += values_[thisCol] * f.dofVector()[ blockNr ][ dofNr ];
          }

          ++ret_it;
        }
      }

      //! return value of entry (row,col)
      field_type operator()(size_type row, size_type col) const
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
      void clearRow(size_type row)
      {
        assert((row>=0) && (row <= dim_[0]));

        nonZeros_[row] = firstCol;
        auto col = row * nz_;
        for(auto i=decltype(nz_){0}; i<nz_; ++i)
        {
          values_ [col] = 0;
          col_[col] = defaultCol;
          ++col;
        }
      }

      //! return max number of non zeros
      //! used in SparseRowMatrixObject::reserve
      size_type numNonZeros() const
      {
        return nz_;
      }

      //! return number of non zeros in row
      //! used in ColCompMatrix::setMatrix
      size_type numNonZeros(size_type i) const
      {
        return nonZeros_[i];
      }

      //! return pair (value,column)
      //! used in ColCompMatrix::setMatrix
      std::pair<const field_type, size_type> realValue(size_type index) const
      {
        return std::pair<const field_type, size_type>(values_[index], col_[index]);
      }

      //! print matrix
      void print(std::ostream& s=std::cout, unsigned int offset=0) const
      {
        std::size_t pos(0);
        for(std::size_t row=0; row<dim_[0]; ++row)
          while(pos < (nz_*(row+1)))
          {
            const auto rv(realValue(pos));
            const auto column(rv.second);
            const auto value(rv.first);
            if((std::abs(value) > 1.e-15) && (column != defaultCol))
              s << row+offset << " " << column+offset << " " << value << std::endl;
            ++pos;
          }
      }

      template <class SizeT, class NumericT >
      void fillCSRStorage( std::vector< std::map<SizeT, NumericT> >& matrix ) const
      {
        matrix.resize( rows() );

        size_type thisCol = 0;
        for(size_type i = 0; i<dim_[ 0 ]; ++i )
        {
          auto& matRow = matrix[ i ];
          for(size_type col=firstCol; col<nz_; ++col, ++thisCol)
          {
            const size_type realCol = col_[thisCol];

            if( realCol == defaultCol )
              continue;

            matRow[ realCol ] = values_[ thisCol ];
          }
        }
      }

    private:
      //! resize matrix
      void resize(size_type rows, size_type cols, size_type nz)
      {
        constexpr auto colVal = defaultCol;
        values_.resize( rows*nz , 0 );
        col_.resize( rows*nz , colVal );
        nonZeros_.resize( rows , 0 );
        dim_[0] = rows;
        dim_[1] = cols;
        nz_ = nz+firstCol;
      }

      //! returns local col index for given global (row,col)
      size_type colIndex(size_type row, size_type col)
      {
        assert((col>=0) && (col <= dim_[1]));
        assert((row>=0) && (row <= dim_[0]));

        size_type i = 0;
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

      std::vector<field_type> values_;
      std::vector<size_type> col_;
      std::vector<size_type> nonZeros_;
      std::array<size_type,2> dim_;
      size_type nz_;
    };



    //! SparseRowMatrixObject
    template< class DomainSpace, class RangeSpace,
              class Matrix = SparseRowMatrix< typename DomainSpace::RangeFieldType > >
    class SparseRowMatrixObject
    {
    protected:
      template< class MatrixObject >
      struct LocalMatrixTraits;

      template< class MatrixObject >
      class LocalMatrix;

    public:
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;
      typedef typename DomainSpaceType::EntityType DomainEntityType;
      typedef typename RangeSpaceType::EntityType RangeEntityType;
      typedef typename DomainSpaceType::EntityType ColumnEntityType;
      typedef typename RangeSpaceType::EntityType RowEntityType;

      typedef typename DomainSpaceType::BlockMapperType DomainBlockMapperType;
      typedef NonBlockMapper< DomainBlockMapperType, DomainSpaceType::localBlockSize > DomainMapperType;
      typedef typename RangeSpaceType::BlockMapperType RangeBlockMapperType;
      typedef NonBlockMapper< RangeBlockMapperType, RangeSpaceType::localBlockSize > RangeMapperType;
      typedef Matrix MatrixType;
      typedef typename MatrixType::size_type size_type;
      typedef typename MatrixType::field_type field_type;

      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, MatrixType > ThisType;

      static const size_type domainLocalBlockSize = DomainSpaceType::dimRange;
      static const size_type rangeLocalBlockSize  = RangeSpaceType::dimRange;

      typedef Dune::FieldMatrix< field_type, rangeLocalBlockSize, domainLocalBlockSize > MatrixBlockType;
      typedef MatrixBlockType  block_type;

      typedef MatrixType PreconditionMatrixType;

      typedef LocalMatrix<ThisType> ObjectType;
      typedef ThisType LocalMatrixFactoryType;
      typedef Fem::ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;
      typedef ColumnObject< ThisType > LocalColumnObjectType;

      //! construct matrix object
      SparseRowMatrixObject( const DomainSpaceType &domainSpace,
                             const RangeSpaceType &rangeSpace,
                             const MatrixParameter& param = SparseRowMatrixParameter() )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        domainMapper_( domainSpace_.blockMapper() ),
        rangeMapper_( rangeSpace_.blockMapper() ),
        sequence_( -1 ),
        matrix_(),
        preconditioning_( param.method() != 0 ),
        localMatrixStack_( *this )
      {}

      //! get domain space (i.e. space that builds the rows)
      const DomainSpaceType& domainSpace() const
      {
        return domainSpace_;
      }

      //! get range space (i.e. space that builds the columns)
      const RangeSpaceType& rangeSpace() const
      {
        return rangeSpace_;
      }

      //! get reference to storage object
      MatrixType &matrix() const
      {
        return matrix_;
      }

      //! interface method from LocalMatrixFactory
      ObjectType *newObject() const
      {
        return new ObjectType( *this, domainSpace_, rangeSpace_, domainMapper_, rangeMapper_ );
      }

      //! get local matrix
      LocalMatrixType localMatrix( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity ) const
      {
        return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
      }

      //! get uninitialized local matrix
      LocalMatrixType localMatrix() const
      {
        return LocalMatrixType( localMatrixStack_ );
      }

      //! get local column
      LocalColumnObjectType localColumn( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType( *this, domainEntity );
      }

      void unitRow( const size_type row )
      {
        for( unsigned int i=0, r = row * domainLocalBlockSize; i<domainLocalBlockSize; ++i, ++r )
        {
          matrix_.clearRow( r );
          matrix_.set( r, r, 1.0 );
        }
      }

      template <class LocalBlock>
      void addBlock( const size_type row, const size_type col, const LocalBlock& block )
      {
        std::array< size_type, rangeLocalBlockSize  > rows;
        std::array< size_type, domainLocalBlockSize > cols;
        for( unsigned int i=0, r = row * domainLocalBlockSize, c = col * domainLocalBlockSize; i<domainLocalBlockSize; ++i, ++r, ++c )
        {
          rows[ i ] = r;
          cols[ i ] = c;
        }

        for( unsigned int i=0; i<domainLocalBlockSize; ++i )
        {
          for( unsigned int j=0; j<domainLocalBlockSize; ++j )
          {
            matrix_.add( rows[ i ], cols[ j ], block[ i ][ j ]);
          }
        }
      }

      template <class LocalBlock>
      void setBlock( const size_type row, const size_type col, const LocalBlock& block )
      {
        std::array< size_type, rangeLocalBlockSize  > rows;
        std::array< size_type, domainLocalBlockSize > cols;
        for( unsigned int i=0, r = row * domainLocalBlockSize, c = col * domainLocalBlockSize; i<domainLocalBlockSize; ++i, ++r, ++c )
        {
          rows[ i ] = r;
          cols[ i ] = c;
        }

        for( unsigned int i=0; i<domainLocalBlockSize; ++i )
        {
          for( unsigned int j=0; j<domainLocalBlockSize; ++j )
          {
            matrix_.set( rows[ i ], cols[ j ], block[ i ][ j ]);
          }
        }
      }

      template< class LocalMatrix >
      void addLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        auto functor = [ &localMat, this ] ( std::pair< int, int > local, const std::pair< size_type, size_type >& global )
        {
          matrix_.add( global.first, global.second, localMat.get( local.first, local.second ) );
        };

        rangeMapper_.mapEach( rangeEntity, makePairFunctor( domainMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat, const Scalar &s )
      {
        auto functor = [ &localMat, &s, this ] ( std::pair< int, int > local, const std::pair< size_type, size_type >& global )
        {
          matrix_.add( global.first, global.second, s * localMat.get( local.first, local.second ) );
        };

        rangeMapper_.mapEach( rangeEntity, makePairFunctor( domainMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix >
      void setLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        auto functor = [ &localMat, this ] ( std::pair< int, int > local, const std::pair< size_type, size_type >& global )
        {
          matrix_.set( global.first, global.second, localMat.get( local.first, local.second ) );
        };

        rangeMapper_.mapEach( rangeEntity, makePairFunctor( domainMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix >
      void getLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, LocalMatrix &localMat ) const
      {
        auto functor = [ &localMat, this ] ( std::pair< int, int > local, const std::pair< size_type, size_type >& global )
        {
          localMat.set( local.first, local.second, matrix_( global.first, global.second ) );
        };

        rangeMapper_.mapEach( rangeEntity, makePairFunctor( domainMapper_, domainEntity, functor ) );
      }

      //! clear matrix
      void clear()
      {
        matrix_.clear();
      }

      void flushAssembly() {}

      template <class Set>
      void reserve (const std::vector< Set >& sparsityPattern )
      {
        reserve( StencilWrapper< DomainSpaceType,RangeSpaceType, Set >( sparsityPattern ) );
      }

      //! reserve memory
      template <class Stencil>
      void reserve(const Stencil &stencil, bool verbose = false )
      {
        if( sequence_ != domainSpace_.sequence() )
        {
          // if empty grid do nothing (can appear in parallel runs)
          if( (domainSpace_.begin() != domainSpace_.end()) && (rangeSpace_.begin() != rangeSpace_.end()) )
          {
            // output some info
            if( verbose )
            {
              std::cout << "Reserve Matrix with (" << rangeSpace_.size() << "," << domainSpace_.size()<< ")" << std::endl;
              std::cout << "Max number of base functions = (" << rangeMapper_.maxNumDofs() << ","
                << domainMapper_.maxNumDofs() << ")" << std::endl;
            }
            // reserve matrix
            const auto nonZeros = std::max( static_cast<size_type>(stencil.maxNonZerosEstimate()*DomainSpaceType::localBlockSize),
                                            matrix_.numNonZeros() );
            matrix_.reserve( rangeSpace_.size(), domainSpace_.size(), nonZeros );
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
      //! this only works for square matrices
      template < class DiscreteFunctionType >
      void extractDiagonal( DiscreteFunctionType& diag ) const
      {
        assert( matrix_.rows() == matrix_.cols() );
        const auto dofEnd = diag.dend();
        size_type row = 0;
        for( auto dofIt = diag.dbegin(); dofIt != dofEnd; ++dofIt, ++row )
        {
          assert( row < matrix_.rows() );
          (*dofIt) = matrix_( row, row );
        }
      }

      template <class Vector>
      void setUnitRows( const Vector &rows )
      {
        const auto &slaveDofs = domainSpace().slaveDofs();
        for (auto r : rows)
        {
          matrix_.clearRow(r);
          matrix_.set(r,r,slaveDofs.isSlave( r )? 0.0 : 1.0);
        }
      }

      //! resort row numbering in matrix to have ascending numbering
      void resort()
      {
        matrix_.resort();
      }

      //! mult method of matrix object used by oem solver
      field_type ddotOEM( const field_type *v, const field_type *w ) const
      {
        typedef AdaptiveDiscreteFunction< DomainSpace > DomainFunctionType;
        DomainFunctionType V( "ddot V", domainSpace_, v );
        DomainFunctionType W( "ddot W", domainSpace_, w );
        return V.scalarProductDofs( W );
      }

      //! mult method of matrix object used by oem solver
      void multOEM( const field_type *arg, field_type *dest ) const
      {
        typedef AdaptiveDiscreteFunction< DomainSpace > DomainFunctionType;
        typedef AdaptiveDiscreteFunction< RangeSpace > RangeFunctionType;

        DomainFunctionType farg( "multOEM arg", domainSpace_, arg );
        RangeFunctionType fdest( "multOEM dest", rangeSpace_, dest );
        apply( farg, fdest );
      }
    protected:
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;
      DomainMapperType domainMapper_ ;
      RangeMapperType rangeMapper_ ;
      int sequence_;
      mutable MatrixType matrix_;
      bool preconditioning_;
      mutable LocalMatrixStackType localMatrixStack_;
    };



    //! LocalMatrixTraits
    template< class DomainSpace, class RangeSpace, class Matrix >
    template< class MatrixObject >
    struct SparseRowMatrixObject< DomainSpace, RangeSpace, Matrix >::LocalMatrixTraits
    {
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      typedef SparseRowMatrixObject< DomainSpaceType, RangeSpaceType, Matrix > SparseRowMatrixObjectType;

      typedef typename SparseRowMatrixObjectType::template LocalMatrix< MatrixObject > LocalMatrixType;

      typedef typename RangeSpaceType::RangeFieldType RangeFieldType;
      typedef RangeFieldType LittleBlockType;

      typedef typename SparseRowMatrixObjectType::DomainMapperType  DomainMapperType;
      typedef typename SparseRowMatrixObjectType::RangeMapperType   RangeMapperType;
    };



    //! LocalMatrix
    template< class DomainSpace, class RangeSpace, class Matrix >
    template< class MatrixObject >
    class SparseRowMatrixObject< DomainSpace, RangeSpace, Matrix >::LocalMatrix
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
      typedef typename MatrixObjectType::MatrixType MatrixType;

      //! type of entries of little blocks
      typedef typename Traits::RangeFieldType RangeFieldType;

      //! type of the DoFs
      typedef RangeFieldType DofType;

      //! type of little blocks
      typedef typename Traits::LittleBlockType LittleBlockType;

      //! type of nonblocked domain mapper
      typedef typename Traits::DomainMapperType DomainMapperType;
      //! type of nonblocked domain mapper
      typedef typename Traits::RangeMapperType RangeMapperType;

      typedef std::vector< typename RangeMapperType::SizeType > RowIndicesType;
      typedef std::vector< typename DomainMapperType::SizeType > ColumnIndicesType;

      //! constructor
      LocalMatrix( const MatrixObjectType &matrixObject,
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
      size_type rows() const
      {
        return rowIndices_.size();
      }

      //! return number of columns
      size_type columns() const
      {
        return columnIndices_.size();
      }

      //! add value to matrix entry
      void add(size_type localRow, size_type localCol, DofType value)
      {
        assert( value == value );
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.add( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! get matrix entry
      DofType get(size_type localRow, size_type localCol) const
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        return matrix_( rowIndices_[ localRow ], columnIndices_[ localCol ] );
      }

      //! set matrix entry to value
      void set(size_type localRow, size_type localCol, DofType value)
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.set( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
      }

      //! set matrix row to zero except diagonla entry
      void unitRow(size_type localRow)
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.unitRow( rowIndices_[ localRow ] );
      }

      //! set matrix row to zero
      void clearRow( size_type localRow )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        matrix_.clearRow( rowIndices_[localRow]);
      }

      //! set matrix column to zero
      void clearCol( size_type localCol )
      {
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.clearCol( columnIndices_[localCol] );
      }

      //! clear all entries belonging to local matrix
      void clear()
      {
        const auto row = rows();
        for(auto i=decltype(row){0}; i < row; ++i )
          matrix_.clearRow( rowIndices_[ i ] );
      }

      //! scale local matrix with a certain value
      void scale( const DofType& value )
      {
        const auto row = rows();
        for(auto i=decltype(row){0}; i < row; ++i )
          matrix_.scaleRow( rowIndices_[ i ] , value );
      }

      //! resort all global rows of matrix to have ascending numbering
      void resort()
      {
        const auto row = rows();
        for(auto i=decltype(row){0}; i < row; ++i )
          matrix_.resortRow( rowIndices_[ i ] );
      }

    protected:
      MatrixType &matrix_;
      const DomainMapperType& domainMapper_;
      const RangeMapperType& rangeMapper_;
      RowIndicesType rowIndices_;
      ColumnIndicesType columnIndices_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPMATRIX_HH
