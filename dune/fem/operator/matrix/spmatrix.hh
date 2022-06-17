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
#include <dune/fem/solver/parameter.hh>

namespace Dune
{

  namespace Fem
  {

    //! SparseRowMatrix
    template <class T, class IndexT = std::size_t,
              class ValuesVector = std::vector< T >,
              class IndicesVector = std::vector< IndexT > >
    class SparseRowMatrix
    {
    public:
      //! matrix field type
      typedef T field_type;
      //! matrix index type
      typedef IndexT size_type;
      typedef SparseRowMatrix<field_type,size_type,ValuesVector,IndicesVector> ThisType;
      //! type of the base matrix
      //! for consistency with ISTLMatrixObject
      typedef ThisType MatrixBaseType;

      static const size_type defaultCol = std::numeric_limits<size_type>::max();
      static const size_type zeroCol    = std::numeric_limits<size_type>::max()-1;
      static const int firstCol = 0;

      SparseRowMatrix(const ThisType& ) = delete;

      //! construct matrix of zero size
      explicit SparseRowMatrix( const bool threading = true ) :
        values_(0), columns_(0), rows_(0), dim_({{0,0}}), maxNzPerRow_(0), compressed_( false ), threading_( threading )
      {}

      //! construct matrix with 'rows' rows and 'cols' columns,
      //! maximum 'nz' non zero values in each row
      SparseRowMatrix(const size_type rows, const size_type cols, const size_type nz, const bool threading = true ) :
        values_(0), columns_(0), rows_(0), dim_({{0,0}}), maxNzPerRow_(0), compressed_( false ), threading_( threading )
      {
        reserve(rows,cols,nz);
      }

      //! reserve memory for given rows, columns and number of non zeros
      void reserve(const size_type rows, const size_type cols, const size_type nz)
      {
        // if( (rows != dim_[0]) || (cols != dim_[1]) || (nz != maxNzPerRow_))
        resize(rows,cols,nz);
        clear();
      }

      //! reserve memory for given rows, columns and number of non zeros
      template <class Stencil>
      void fillPattern(const Stencil& stencil,
                       const size_type rowBlockSize,
                       const size_type colBlockSize )
      {
        const auto& sparsityPattern = stencil.globalStencil();
        for( const auto& entry : sparsityPattern )
        {
          const auto& blockRow = entry.first;
          const auto& blockColumnSet = entry.second;

          // blocking of rows
          const size_type nextRow = (blockRow + 1) * rowBlockSize;
          for( size_type row = blockRow * rowBlockSize; row < nextRow; ++row )
          {
            size_type column = startRow( row );
            for( const auto& blockColEntry : blockColumnSet )
            {
              size_type col = blockColEntry * colBlockSize;
              for( size_type c = 0; c<colBlockSize; ++c, ++col, ++column )
              {
                assert( column < endRow( row ) );
                columns_[ column ] = col ;
              }
            }
          }
        }
      }

      //! return number of rows
      size_type rows () const
      {
        return dim_[0];
      }

      //! return number of columns
      size_type cols () const
      {
        return dim_[1];
      }

      //! set entry to value (also setting 0 will result in an entry)
      void set(const size_type row, const size_type col, const field_type val)
      {
        assert((col>=0) && (col < dim_[1]));
        assert((row>=0) && (row < dim_[0]));

        const size_type column = colIndex(row,col) ;
        assert( column != defaultCol && column != zeroCol );

        values_ [ column ] = val;
        columns_[ column ] = col;
      }

      //! add value to row,col entry
      void add(const size_type row, const size_type col, const field_type val)
      {
        assert((col>=0) && (col < dim_[1]));
        assert((row>=0) && (row < dim_[0]));

        const size_type column = colIndex(row,col) ;
        assert( column != defaultCol && column != zeroCol );

        values_ [ column ] += val;
        columns_[ column ] = col;
      }

      //! ret = A*f
      template<class ArgDFType, class DestDFType>
      void apply(const ArgDFType& f, DestDFType& ret ) const
      {
        bool runParallel = threading_;

        auto doApply = [this, &f, &ret]()
        {
          constexpr auto blockSize = ArgDFType::DiscreteFunctionSpaceType::localBlockSize;

          // compute slice of rows to be worked on
          const auto slice = sliceBeginEnd( MPIManager::thread(), MPIManager::numThreads(), std::true_type() );

          // same as begin just with a row not necessarily zero
          auto ret_it = ret.dofVector().find( slice.first );

          const auto& fDofVec = f.dofVector();
          for(size_type row = slice.first; row<slice.second; ++row)
          {
            const size_type endrow = endRow( row );
            (*ret_it) = 0.0;
            for(size_type col = startRow( row ); col<endrow; ++col)
            {
              const auto realCol = columns_[ col ];

              if( ! compressed_ && ((realCol == defaultCol) || (realCol == zeroCol)) )
                continue;

              const auto blockNr = realCol / blockSize ;
              const auto dofNr   = realCol % blockSize ;
              (*ret_it) += values_[ col ] * fDofVec[ blockNr ][ dofNr ];
            }

            ++ret_it;
          }
        };

        if( runParallel )
        {
          try {
            // execute in parallel
            MPIManager :: run ( doApply );
          }
          catch ( const SingleThreadModeError& e )
          {
            runParallel = false;
          }
        }

        // run serial if threading disabled or something else went wrong
        if( ! runParallel )
        {
          doApply();
        }
      }

      //! return value of entry (row,col)
      field_type get(const size_type row, const size_type col) const
      {
        assert((col>=0) && (col < dim_[1]));
        assert((row>=0) && (row < dim_[0]));

        const size_type endrow = endRow( row );
        for( size_type i = startRow( row ); i<endrow; ++i )
        {
          if(columns_[ i ] == col)
            return values_[ i ];
        }
        return 0;
      }

      //! set all matrix entries to zero
      void clear()
      {
        std::fill( values_.begin(), values_.end(), 0 );
        for (auto &c : columns_) c = defaultCol;
      }

      //! set all entries in row to zero
      void clearRow(const size_type row)
      {
        assert((row>=0) && (row < dim_[0]));

        const size_type endrow = endRow( row );
        for(size_type idx = startRow( row ); idx<endrow; ++idx )
        {
          values_[ idx ] = 0;
          // if ( !compressed_ )
          //    columns_[idx] = zeroCol;
        }
      }

      //! scale all entries in row with a given value
      void scale(const size_type row, const size_type col, const field_type val)
      {
        assert((row>=0) && (row < rows()) );
        assert((col>=0) && (col < cols()) );

        const size_type column = colIndex(row,col) ;
        assert( column != defaultCol && column != zeroCol );

        // scale value
        values_ [ column ] *= val;
      }

      //! return max number of non zeros
      //! used in SparseRowMatrixObject::reserve
      size_type maxNzPerRow() const
      {
        return maxNzPerRow_;
      }

      //! return max number of non zeros
      //! used in SparseRowMatrixObject::reserve
      size_type numNonZeros() const
      {
        return dim_[0] > 0 ? rows_[ dim_[0] ] : 0;
      }

      //! return number of non zeros in row
      //! used in ColCompMatrix::setMatrix
      size_type numNonZeros(size_type row) const
      {
        assert( row >= 0 && row < dim_[0] );
        return endRow( row ) - startRow( row );
      }

      //! return pair (value,column)
      //! used in ColCompMatrix::setMatrix and FemPy CRS matrix export
      std::pair<const field_type, size_type> realValue(size_type index) const
      {
        return std::pair<const field_type, size_type>(values_[index], columns_[index]);
      }

      //! print matrix
      void print(std::ostream& s=std::cout, unsigned int offset=0) const
      {
        for(size_type row=0; row<dim_[0]; ++row)
        {
          const size_type endrow = endRow( row );
          for( size_type pos = startRow( row ); pos<endrow; ++pos )
          {
            const auto rv(realValue(pos));
            const auto column(rv.second);
            const auto value(rv.first);
            if((std::abs(value) > 1.e-15) && (column != defaultCol))
              s << row+offset << " " << column+offset << " " << value << std::endl;
          }
        }
      }

      template <class SizeT, class NumericT >
      void fillCSRStorage( std::vector< std::map<SizeT, NumericT> >& matrix ) const
      {
        matrix.resize( dim_[0] );

        size_type thisCol = 0;
        for(size_type row = 0; row<dim_[ 0 ]; ++row )
        {
          auto& matRow = matrix[ row ];
          const size_type endrow = endRow( row );
          for(size_type col = startRow( row ); col<endrow; ++col)
          {
            const size_type realCol = columns_[ col ];

            if( ! compressed_ && (realCol == defaultCol || realCol == zeroCol) )
              continue;

            matRow[ realCol ] = values_[ thisCol ];
          }
        }
      }

      void compress ()
      {
        if( ! compressed_ && (dim_[0] != 0) && (dim_[1] != 0))
        {
          // determine first row nonZeros
          size_type newpos = 0 ;
          for( newpos = startRow( 0 ); newpos < endRow( 0 ); ++newpos )
          {
            if( columns_[ newpos ] == defaultCol )
              break;
          }

          for( size_type row = 1; row<dim_[0]; ++row )
          {
            const size_type endrow = endRow( row );
            size_type col = startRow( row );
            // update new row start position
            rows_[ row ] = newpos;
            for(; col<endrow; ++col, ++newpos )
            {
              if( columns_[ col ] == defaultCol )
                break ;

              assert( newpos <= col );
              values_ [ newpos ] = values_ [ col ];
              columns_[ newpos ] = columns_[ col ];
            }
          }
          rows_[ dim_[0] ] = newpos ;

          // values_.resize( newpos );
          // columns_.resize( newpos );
          compressed_ = true ;
        }
      }

      size_type startRow ( const size_type row ) const
      {
        return rows_[ row ];
      }

      size_type endRow ( const size_type row ) const
      {
        return rows_[ row+1 ];
      }

      std::tuple< ValuesVector&, IndicesVector&, IndicesVector& > exportCRS()
      {
        // only return internal data in compressed status
        compress();
        return std::tie(values_,columns_,rows_);
      }

      //! Apply Jacobi/SOR method
      template<class DiagType, class ArgDFType, class DestDFType, class WType>
      void forwardIterative(const DiagType& diagInv, const ArgDFType& b, const DestDFType& xold, DestDFType& xnew, const WType& w ) const
      {
        parallelIterative(diagInv, b, xold, xnew, w, std::true_type() );
      }

      //! Apply Jacobi/SOR method
      template<class DiagType, class ArgDFType, class DestDFType, class WType>
      void backwardIterative(const DiagType& diagInv, const ArgDFType& b, const DestDFType& xold, DestDFType& xnew, const WType& w ) const
      {
        parallelIterative(diagInv, b, xold, xnew, w, std::false_type() );
      }

    protected:
      // return slice of rows for given thread in forward direction
      std::pair< size_type, size_type > sliceBeginEnd(const size_type thread,  const size_type numThreads, std::true_type ) const
      {
        const size_type sliceSize  = this->rows() / numThreads;
        const size_type sliceBegin = (thread * sliceSize) ;
        const size_type sliceEnd   = (thread == numThreads-1 ) ? this->rows(): (sliceBegin + sliceSize);
        return std::make_pair( sliceBegin, sliceEnd );
      }

      // return slice of rows for given thread in backward direction
      std::pair< size_type, size_type > sliceBeginEnd(const size_type thread, const size_type numThreads, std::false_type ) const
      {
        std::pair< size_type, size_type > beginEnd = sliceBeginEnd( thread, numThreads, std::true_type() );
        return std::make_pair( beginEnd.second-1, beginEnd.first-1 );
      }

      template<class DiagType, class ArgDFType, class DestDFType, class WType, bool forward >
      void parallelIterative(const DiagType& diagInv, const ArgDFType& b, const DestDFType& xold, DestDFType& xnew,
                             const WType& w, std::integral_constant<bool, forward> direction ) const
      {
        bool runParallel = threading_ && (&xold != &xnew) ; // only for Jacobi

        auto doIterate = [this, &diagInv, &b, &xold, &xnew, &w, &direction] ()
        {
          // compute slice to be worked on depending on direction
          const auto slice = sliceBeginEnd( MPIManager::thread(), MPIManager::numThreads(), direction );

          doParallelIterative( diagInv.dofVector().find( slice.first ), // still O(1) operation just like for begin()
                               b.dofVector().find( slice.first ),
                               xold,
                               xnew.dofVector().find( slice.first ),
                               w,
                               slice.first, // row begin
                               slice.second,   // row end
                               direction );
        };

        if( runParallel )
        {
          try {
            // execute in parallel
            MPIManager :: run ( doIterate );
          }
          catch ( const SingleThreadModeError& e )
          {
            runParallel = false;
          }
        }

        // run serial if threading disabled or something else went wrong
        if( ! runParallel )
        {
          doIterate();
        }
      }

      //! Apply Jacobi/SOR method
      template<class DiagIt, class ArgDFIt, class DestDFType, class DestDFIt,
               class WType, bool forward>
      void doParallelIterative(DiagIt diag, ArgDFIt bit, const DestDFType& xold, DestDFIt xit,
                               const WType& w,
                               size_type row, const size_type end,
                               std::integral_constant<bool, forward> ) const
      {
        constexpr auto blockSize = DestDFType::DiscreteFunctionSpaceType::localBlockSize;

        const auto nextRow = [&diag, &xit, &bit](size_type &row)
        {
          if constexpr ( forward )
          {
            ++diag; ++xit; ++bit; ++row;
          }
          else
          {
            --diag; --xit; --bit; --row;
          }
        };

        const auto& xOldVec = xold.dofVector();
        for(; row != end; nextRow(row))
        {
          auto rhs = (*bit);

          const size_type endrow = endRow( row );
          for(size_type col = startRow( row ); col<endrow; ++col)
          {
            const auto realCol = columns_[ col ];

            if( (realCol == row ) || (! compressed_ && ((realCol == defaultCol) || (realCol == zeroCol))) )
              continue;

            const auto blockNr = realCol / blockSize ;
            const auto dofNr   = realCol % blockSize ;

            rhs -= values_[ col ] * xOldVec[ blockNr ][ dofNr ] ;
          }

          (*xit) = w * (rhs * (*diag));
        }
      }
      //! resize matrix
      void resize(size_type rows, size_type cols, size_type nz)
      {
        constexpr auto colVal = defaultCol;
        values_.resize( rows*nz , 0 );
        columns_.resize( rows*nz , colVal );
        rows_.resize( rows+1 , 0 );
        rows_[ 0 ] = 0;
        for( size_type i=1; i <= rows; ++i )
        {
          rows_[ i ] = rows_[ i-1 ] + nz ;
        }
        compressed_ = false;

        dim_[0] = rows;
        dim_[1] = cols;
        maxNzPerRow_ = nz+firstCol;
      }

      //! returns local col index for given global (row,col)
      size_type colIndex(size_type row, size_type col)
      {
        assert((col>=0) && (col < dim_[1]));
        assert((row>=0) && (row < dim_[0]));

        const size_type endR  = endRow( row );
        size_type i = startRow( row );
        // find local column or empty spot
        for( ;  i < endR; ++i )
        {
          if( columns_[ i ] == defaultCol || columns_[ i ] == zeroCol || columns_[ i ] == col )
          {
            return i;
          }
        }

        assert(0);
        DUNE_THROW( InvalidStateException, "Could not store entry in sparse matrix - no space available" );

        // TODO: implement resize with 2*nz
        std::abort();
        return defaultCol;

        /*
        if(columns_[ i ] == col)
          return i;  // column already in matrix
        else if( columns_[ i ] == defaultCol )
        { // add this column at end of this row
          ++nonZeros_[row];
          return i;
        }
        else
        {
          std::abort();
          // TODO re-implement
          //
          ++nonZeros_[row];
          // must shift this row to add col at the position i
          auto j = nz_-1; // last column
          if (columns_[row*nz_+j] != defaultCol)
          { // new space available - so resize
            resize( rows(), cols(), (2 * nz_) );
            j++;
          }
          for(;j>i;--j)
          {
            columns_[row*nz_+j] = columns_[row*nz_+j-1];
            values_[row*nz_+j] = values_[row*nz_+j-1];
          }
          columns_[row*nz_+i] = col;
          values_[row*nz_+i] = 0;
          return i;
        }
        */
      }

      ValuesVector  values_;
      IndicesVector columns_;
      IndicesVector rows_;

      std::array<size_type,2> dim_;
      size_type maxNzPerRow_;
      bool compressed_;
      const bool threading_;
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
                             const SolverParameter& param = SolverParameter() )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        domainMapper_( domainSpace_.blockMapper() ),
        rangeMapper_( rangeSpace_.blockMapper() ),
        sequence_( -1 ),
        matrix_( param.threading() ),
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

    protected:
      //! get reference to storage object, for internal use
      MatrixType& matrix() const
      {
        return matrix_;
      }

      void finalizeAssembly() const { const_cast< ThisType& > (*this).compress(); }

    public:
      //! get reference to storage object
      MatrixType& exportMatrix() const
      {
        finalizeAssembly();
        return matrix_;
      }


      //! interface method from LocalMatrixFactory
      ObjectType* newObject() const
      {
        return new ObjectType( *this, domainSpace_, rangeSpace_, domainMapper_, rangeMapper_ );
      }

      /** \deprecated Use TemporaryLocalMatrix in combination with
        *             {add,set,get}LocalMatrix on matrix object
        *  return local matrix object
        */
      [[deprecated("Use TemporaryLocalMatrix,LocalContribution and {get,add,set}LocalMatrix")]]
      LocalMatrixType localMatrix( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity ) const
      {
        return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
      }

      /** \deprecated Use TemporaryLocalMatrix in combination with
        *             {add,set,get}LocalMatrix on matrix object
        *  return local matrix object
        */
      [[deprecated("Use TemporaryLocalMatrix,LocalContribution and {get,add,set}LocalMatrix")]]
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
          localMat.set( local.first, local.second, matrix_.get( global.first, global.second ) );
        };

        rangeMapper_.mapEach( rangeEntity, makePairFunctor( domainMapper_, domainEntity, functor ) );
      }

      //! clear matrix
      void clear()
      {
        matrix_.clear();
      }

      //! compress matrix to a real CRS format
      void compress() { matrix_.compress(); }

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
                                            matrix_.maxNzPerRow() );
            matrix_.reserve( rangeSpace_.size(), domainSpace_.size(), nonZeros );
            matrix_.fillPattern( stencil, RangeSpaceType::localBlockSize, DomainSpaceType::localBlockSize );
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
          (*dofIt) = matrix_.get( row, row );
        }
      }

      template <class Container>
      void setUnitRows( const Container& unitRows, const Container& auxRows )
      {
        for (const auto& r : unitRows )
        {
          matrix_.clearRow(r);
          matrix_.set(r, r, 1.0);
        }

        for (const auto& r : auxRows )
        {
          matrix_.clearRow(r);
          // not sure if this is really needed,
          // but for consistency with previous code
          matrix_.set(r, r, 0.0);
        }
      }

      //! resort row numbering in matrix to have ascending numbering
      void resort()
      {
        DUNE_THROW(NotImplemented,"SpMatrixObject::resort is not implemented");
        // this method does not even exist on SpMatrix!!!
        // matrix_.resort();
      }

    protected:
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;
      DomainMapperType domainMapper_ ;
      RangeMapperType rangeMapper_ ;
      int sequence_;
      mutable MatrixType matrix_;
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
        // columns are determined by the domain space
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
        return matrix_.get( rowIndices_[ localRow ], columnIndices_[ localCol ] );
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
        const size_type nrows = rows();
        for(size_type i=0; i < nrows; ++i )
          matrix_.clearRow( rowIndices_[ i ] );
      }

      //! resort all global rows of matrix to have ascending numbering
      void resort()
      {
        DUNE_THROW(NotImplemented,"SpMatrixObject::LocalMatrix::resort is not implemented");
        //const size_type nrows = rows();
        //for(size_type i=0; i < nrows; ++i )
          //matrix_.resortRow( rowIndices_[ i ] );
      }

      //! scale local matrix with a certain value
      void scale( const DofType& value )
      {
        const size_type nrows = rows();
        const size_type ncols = columns();
        for(size_type i=0; i < nrows; ++i )
        {
          for( size_type j=0; j < ncols; ++j )
          {
            scale(i, j, value );
          }
        }
      }

    protected:
      //! scale matrix entry with value
      void scale(size_type localRow, size_type localCol, DofType value)
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        matrix_.scale( rowIndices_[ localRow ], columnIndices_[ localCol ], value );
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
