#ifndef DUNE_FEM_ISTLMATRIXWRAPPER_HH
#define DUNE_FEM_ISTLMATRIXWRAPPER_HH

#if HAVE_DUNE_ISTL

//- system includes
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <string>

//- Dune common includes
#include <dune/common/exceptions.hh>
#include <dune/common/fmatrix.hh>

//- Dune istl includes
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

//- Dune fem includes
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/storage/objectstack.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/operator/matrix/istlpreconditioner.hh>
#include <dune/fem/operator/matrix/functor.hh>

namespace Dune
{

  namespace Fem
  {


    //! forward declarations
    template< class MatrixObject >
    class ISTLLocalMatrix;

    template <class DomainSpaceImp, class RangeSpaceImp,
              class DomainBlock = Dune::FieldVector< typename DomainSpaceImp :: RangeFieldType, DomainSpaceImp:: localBlockSize >,
              class RangeBlock  = Dune::FieldVector< typename RangeSpaceImp  :: RangeFieldType, RangeSpaceImp :: localBlockSize > >
    class ISTLMatrixObject;

    ///////////////////////////////////////////////////////
    // --ISTLMatrixHandle
    //////////////////////////////////////////////////////
    template <class LittleBlockType, class RowDiscreteFunctionImp, class ColDiscreteFunctionImp = RowDiscreteFunctionImp>
    class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType>
    {
        friend struct MatrixDimension<ImprovedBCRSMatrix>;
      public:
        typedef RowDiscreteFunctionImp RowDiscreteFunctionType;
        typedef ColDiscreteFunctionImp ColDiscreteFunctionType;

        typedef BCRSMatrix<LittleBlockType> BaseType;
        //! type of the base matrix
        typedef BaseType MatrixBaseType;
        typedef typename BaseType :: RowIterator RowIteratorType ;
        typedef typename BaseType :: ColIterator ColIteratorType ;

        typedef ImprovedBCRSMatrix< LittleBlockType, RowDiscreteFunctionImp, ColDiscreteFunctionImp > ThisType;

        typedef typename BaseType :: size_type size_type;

        //===== type definitions and constants

        //! export the type representing the field
        typedef typename BaseType::field_type field_type;

        //! export the type representing the components
        typedef typename BaseType::block_type block_type;

        //! export the allocator type
        typedef typename BaseType:: allocator_type allocator_type;

        //! implement row_type with compressed vector
        typedef typename BaseType :: row_type row_type;

        /** \brief Iterator for the entries of each row */
        typedef typename BaseType :: ColIterator ColIterator;

        /** \brief Iterator for the entries of each row */
        typedef typename BaseType :: ConstColIterator ConstColIterator;

        /** \brief Const iterator over the matrix rows */
        typedef typename BaseType :: RowIterator RowIterator;

        /** \brief Const iterator over the matrix rows */
        typedef typename BaseType :: ConstRowIterator ConstRowIterator;

        //! type of discrete function space
        typedef typename ColDiscreteFunctionType :: DiscreteFunctionSpaceType RangeSpaceType;

        //! type of row block vector
        typedef typename RowDiscreteFunctionType :: DofStorageType  RowBlockVectorType;

        //! type of column block vector
        typedef typename ColDiscreteFunctionType :: DofStorageType  ColBlockVectorType;

        //! type of communication object
        typedef typename RangeSpaceType :: GridType :: Traits :: Communication   CommunicationType ;

        typedef typename BaseType :: BuildMode BuildMode ;

      public:
        //! constructor used by ISTLMatrixObject to build matrix in implicit mode
        ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nnz, double overflowFraction) :
          BaseType (rows, cols, nnz, overflowFraction, BaseType::implicit)
        {}

        //! constructor using old row_wise assembly (and used by ILU preconditioner)
        ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nz = 0 ) :
          BaseType (rows, cols, BaseType :: row_wise)
        {}

        //! copy constructor, needed by ISTL preconditioners
        ImprovedBCRSMatrix( ) :
          BaseType ()
        {}

        //! copy constructor, needed by ISTL preconditioners
        ImprovedBCRSMatrix(const ImprovedBCRSMatrix& org) :
          BaseType(org)
        {}

        ConstRowIterator slicedBegin( const size_type row ) const
        {
          return ConstRowIterator( this->r, row );
        }

        ConstRowIterator slicedEnd( const size_type row ) const
        {
          return ConstRowIterator( this->r, row );
        }

        std::pair< size_type, size_type > sliceBeginEnd( const size_type thread, const size_type numThreads ) const
        {
          const size_type sliceSize  = this->N() / numThreads ;
          const size_type sliceBegin = thread * sliceSize ;
          const size_type sliceEnd   = ( thread == numThreads-1 ) ? this->N() : (sliceBegin + sliceSize) ;
          return std::make_pair( sliceBegin, sliceEnd );
        }

      protected:
        template <class X, class Y, bool isMV>
        void mvThreadedImpl( const field_type& alpha,
                             const X& x, Y& y, std::integral_constant<bool, isMV> ) const
        {
          auto doMV = [this, &alpha, &x, &y] ()
          {
            const auto slice = sliceBeginEnd( MPIManager::thread(), MPIManager::numThreads() );
            const ConstRowIterator endi = slicedEnd( slice.second );
            for (ConstRowIterator i = slicedBegin(slice.first); i!=endi; ++i)
            {
              if constexpr ( isMV )
                y[i.index()] = 0;

              ConstColIterator endj = (*i).end();
              for (ConstColIterator j=(*i).begin(); j!=endj; ++j)
              {
                auto&& xj = Dune::Impl::asVector(x[j.index()]);
                auto&& yi = Dune::Impl::asVector(y[i.index()]);
                if constexpr ( isMV )
                  Dune::Impl::asMatrix(*j).umv(xj, yi);
                else
                  Dune::Impl::asMatrix(*j).usmv(alpha, xj, yi);
              }
            }
          };

          try {
            // execute in parallel
            MPIManager :: run ( doMV );
          }
          catch ( const SingleThreadModeError& e )
          {
            // serial version
            if constexpr ( isMV )
              this->mv( x, y );
            else
            {
              DUNE_THROW(SingleThreadModeError,"ImprovedBCRSMatrix::usmvThreaded: Cannot recover from threading error. Disable threading!");
              // this->usmv( alpha, x, y );
            }
          }
        }

      public:
        template <class X, class Y>
        void usmvThreaded( const field_type& alpha, const X& x, Y& y ) const
        {
          mvThreadedImpl(alpha, x, y, std::false_type() );
        }

        template <class X, class Y>
        void mvThreaded( const X& x, Y& y ) const
        {
          mvThreadedImpl(1.0, x, y, std::true_type() );
        }

        template < class SparsityPattern >
        void createEntries(const SparsityPattern& sparsityPattern)
        {
          const auto spEnd = sparsityPattern.end();

          // insert map of indices into matrix
          auto endcreate = this->createend();
          for(auto create = this->createbegin(); create != endcreate; ++create)
          {
            const auto row = sparsityPattern.find( create.index() );
            // if a row is empty then do nothing
            if( row == spEnd ) continue;

            // insert all indices for this row
            for ( const auto& col : row->second )
              create.insert( col );
          }
        }

        //! clear Matrix, i.e. set all entries to 0
        void clear()
        {
          for (auto& row : *this)
            for (auto& entry : row)
              entry = 0;
        }

        //! clear Matrix, i.e. set all entries to 0
        void unitRow( const size_t row )
        {
          block_type idBlock( 0 );
          for (int i = 0; i < idBlock.rows; ++i)
              idBlock[i][i] = 1.0;

          auto& matRow = (*this)[ row ];
          auto colIt = matRow.begin();
          const auto& colEndIt = matRow.end();
          for (; colIt != colEndIt; ++colIt)
          {
              if( colIt.index() == row )
                  *colIt = idBlock;
              else
                  *colIt = 0.0;
          }
        }

        //! setup like the old matrix but remove rows with hanging nodes
        template <class HangingNodesType>
        void setup(ThisType& oldMatrix, const HangingNodesType& hangingNodes)
        {
          // necessary because element traversal not necessarily is in ascending order
          typedef std::set< std::pair<int, block_type> > LocalEntryType;
          typedef std::map< int , LocalEntryType > EntriesType;
          EntriesType entries;

          // map of indices
          std::map< int , std::set<int> > indices;
          // not insert map of indices into matrix
          auto rowend  = oldMatrix.end();
          for(auto it  = oldMatrix.begin(); it != rowend; ++it)
          {
            const auto row = it.index();
            auto& localIndices = indices[ row ];

            if( hangingNodes.isHangingNode( row ) )
            {
              // insert columns into other columns
              const auto& cols = hangingNodes.associatedDofs( row );
              const auto colSize = cols.size();
              for(auto i=0; i<colSize; ++i)
              {
                assert( ! hangingNodes.isHangingNode( cols[i].first ) );

                // get local indices of col
                auto& localColIndices = indices[ cols[i].first ];
                auto& localEntry = entries[  cols[i].first ];

                // copy from old matrix
                auto endj = (*it).end();
                for (auto j= (*it).begin(); j!=endj; ++j)
                {
                  localColIndices.insert( j.index () );
                  localEntry.insert( std::make_pair( j.index(), (cols[i].second * (*j)) ));
                }
              }

              // insert diagonal and hanging columns
              localIndices.insert( row );
              for(auto i=0; i<colSize; ++i)
                localIndices.insert( cols[i].first );
            }
            else
            {
              // copy from old matrix
              auto endj = (*it).end();
              for (auto j= (*it).begin(); j!=endj; ++j)
                localIndices.insert( j.index () );
            }
          }

          // create matrix from entry map
          createEntries( indices );

          // not insert map of indices into matrix
          auto rowit  = oldMatrix.begin();

          auto endcreate = this->end();
          for(auto create = this->begin(); create != endcreate; ++create, ++rowit )
          {
            assert( rowit != oldMatrix.end() );

            const auto row = create.index();
            if( hangingNodes.isHangingNode( row ) )
            {
              const auto& cols = hangingNodes.associatedDofs( row );

              std::map< const int , block_type > colMap;
              // only working for block size 1 ath the moment
              assert( block_type :: rows == 1 );
              // insert columns into map
              const auto colSize = cols.size();
              for( auto i=0; i<colSize; ++i)
                colMap[ cols[i].first ] = -cols[i].second;
              // insert diagonal into map
              colMap[ row ] = 1;

              auto endj = (*create).end();
              for (auto j= (*create).begin(); j!=endj; ++j)
              {
                assert( colMap.find( j.index() ) != colMap.end() );
                (*j) = colMap[ j.index() ];
              }
            }
            // if entries are equal, just copy
            else if ( entries.find( row ) == entries.end() )
            {
              auto colit = (*rowit).begin();
              auto endj = (*create).end();
              for (auto j= (*create).begin(); j!=endj; ++j, ++colit )
              {
                assert( colit != (*rowit).end() );
                (*j) = (*colit);
              }
            }
            else
            {
              std::map< int , block_type > oldCols;

              {
                auto colend = (*rowit).end();
                for(auto colit = (*rowit).begin(); colit != colend; ++colit)
                  oldCols[ colit.index() ] = 0;
              }

              auto entry = entries.find( row );
              assert( entry  != entries.end ());

              {
                auto endcol = (*entry).second.end();
                for( auto co = (*entry).second.begin(); co != endcol; ++co)
                  oldCols[ (*co).first ] = 0;
              }

              {
                auto colend = (*rowit).end();
                for(auto colit = (*rowit).begin(); colit != colend; ++colit)
                  oldCols[ colit.index() ] += (*colit);
              }

              {
                auto endcol = (*entry).second.end();
                for( auto co = (*entry).second.begin(); co != endcol; ++co)
                  oldCols[ (*co).first ] += (*co).second;
              }

              auto endj = (*create).end();
              for (auto j= (*create).begin(); j!=endj; ++j )
              {
                auto colEntry = oldCols.find( j.index() );
                if( colEntry != oldCols.end() )
                  (*j) = (*colEntry).second;
                else
                  abort();
              }
            }
          }
        }

        //! extract diagonal of matrix to block vector
        void extractDiagonal( ColBlockVectorType& diag ) const
        {
          const auto endi = this->end();
          for (auto i = this->begin(); i!=endi; ++i)
          {
            // get diagonal entry of matrix
            const auto row = i.index();
            auto entry = (*i).find( row );
            const LittleBlockType& block = (*entry);
            enum { blockSize = LittleBlockType :: rows };
            for( auto l=0; l<blockSize; ++l )
              diag[ row ][ l ] = block[ l ][ l ];
          }
        }

        //! return value of entry (row,col) where row and col are global indices not block wise
        //! in order to be consistent with SparseRowMatrix.
        field_type operator()(const std::size_t row, const std::size_t col) const
        {
          const std::size_t blockRow(row/(LittleBlockType :: rows));
          const std::size_t localRowIdx(row%(LittleBlockType :: rows));
          const std::size_t blockCol(col/(LittleBlockType :: cols));
          const std::size_t localColIdx(col%(LittleBlockType :: cols));

          const auto& matrixRow(this->operator[](blockRow));
          auto entry = matrixRow.find( blockCol );
          const LittleBlockType& block = (*entry);
          return block[localRowIdx][localColIdx];
        }

        //! set entry to value (row,col) where row and col are global indices not block wise
        //! in order to be consistent with SparseRowMatrix.
        void set(const std::size_t row, const std::size_t col, field_type value)
        {
          const std::size_t blockRow(row/(LittleBlockType :: rows));
          const std::size_t localRowIdx(row%(LittleBlockType :: rows));
          const std::size_t blockCol(col/(LittleBlockType :: cols));
          const std::size_t localColIdx(col%(LittleBlockType :: cols));

          auto& matrixRow(this->operator[](blockRow));
          auto entry = matrixRow.find( blockCol );
          LittleBlockType& block = (*entry);
          block[localRowIdx][localColIdx] = value;
        }

        //! print matrix
        void print(std::ostream& s=std::cout, unsigned int offset=0) const
        {
          s.precision( 6 );
          const auto endi=this->end();
          for (auto i=this->begin(); i!=endi; ++i)
          {
            const auto endj = (*i).end();
            for (auto j=(*i).begin(); j!=endj; ++j)
              if( (*j).infinity_norm() > 1.e-15)
                s << i.index()+offset << " " << j.index()+offset << " " << *j << std::endl;
          }
        }
    };



    // ISTLLocalMatrixTraits
    // ---------------------

    template< class MatrixObject >
    struct ISTLLocalMatrixTraits
    {
      typedef typename MatrixObject::DomainSpaceType DomainSpaceType;
      typedef typename MatrixObject::RangeSpaceType RangeSpaceType;
      typedef typename DomainSpaceType::RangeFieldType RangeFieldType;

      typedef ISTLLocalMatrix< MatrixObject > LocalMatrixType;
      typedef typename MatrixObject::MatrixType::block_type LittleBlockType;
    };


    // ISTLLocalMatrix
    // ---------------

    template< class MatrixObject >
    class ISTLLocalMatrix
    : public LocalMatrixDefault< ISTLLocalMatrixTraits< MatrixObject > >
    {
    public:
      //! type of base class
      typedef LocalMatrixDefault< ISTLLocalMatrixTraits< MatrixObject > > BaseType;

      //! type of matrix object
      typedef MatrixObject MatrixObjectType;
      //! type of matrix
      typedef typename MatrixObjectType::MatrixType MatrixType;
      //! type of little blocks
      typedef typename MatrixType::block_type LittleBlockType;
      // type of row and column indices
      typedef typename MatrixType :: size_type size_type;

      typedef typename MatrixObjectType::DomainSpaceType DomainSpaceType;
      typedef typename MatrixObjectType::RangeSpaceType RangeSpaceType;

      typedef typename MatrixObjectType::DomainEntityType DomainEntityType;
      typedef typename MatrixObjectType::RangeEntityType RangeEntityType;

      //! type of entries of little blocks
      typedef typename DomainSpaceType::RangeFieldType DofType;
      typedef typename MatrixType::row_type RowType;

      //! type of row mapper
      typedef typename DomainSpaceType::BlockMapperType ColMapperType;
      //! type of col mapper
      typedef typename RangeSpaceType::BlockMapperType RowMapperType;

      static const int littleCols = MatrixObjectType::littleCols;
      static const int littleRows = MatrixObjectType::littleRows;

      typedef typename MatrixType::size_type Index;

      /** \deprecated Use TemporaryLocalMatrix in combination with
        *             {add,set,get}LocalMatrix on matrix object
        */
      [[deprecated("Use TemporaryLocalMatrix,LocalContribution and {get,add,set}LocalMatrix")]]
      ISTLLocalMatrix ( const MatrixObjectType& mObj, const DomainSpaceType& domainSpace, const RangeSpaceType& rangeSpace )
        : BaseType( domainSpace, rangeSpace ),
          rowMapper_( rangeSpace.blockMapper() ),
          colMapper_( domainSpace.blockMapper() ),
          numRows_( rowMapper_.maxNumDofs() ),
          numCols_( colMapper_.maxNumDofs() ),
          matrixObj_( mObj )
      {}

      /** \deprecated Use TemporaryLocalMatrix in combination with
        *             {add,set,get}LocalMatrix on matrix object
        */
      [[deprecated("Use TemporaryLocalMatrix,LocalContribution and {get,add,set}LocalMatrix")]]
      ISTLLocalMatrix ( const ISTLLocalMatrix& org )
        : BaseType( org ),
          rowMapper_(org.rowMapper_),
          colMapper_(org.colMapper_),
          numRows_( org.numRows_ ),
          numCols_( org.numCols_ ),
          matrices_(org.matrices_),
          matrixObj_(org.matrixObj_)
      {}

      //! initialize this local Matrix to (colEntity, rowEntity)
      void init ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        // initialize base functions sets
        BaseType :: init ( domainEntity, rangeEntity );

        numRows_  = rowMapper_.numDofs( rangeEntity );
        numCols_  = colMapper_.numDofs( domainEntity );

        // resize matrix pointer storage for new row/col numbers
        matrices_.resize( numRows_ );
        for( auto& row : matrices_ )
        {
          row.resize( numCols_, nullptr );
        }

        if(  matrixObj_.implicitModeActive() )
        {
          auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
          {
            return matrixObj_.matrix().entry( index.first, index.second );
          };
          initBlocks( blockAccess, domainEntity, rangeEntity );
        }
        else
        {
          auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
          {
            return matrixObj_.matrix()[ index.first ][ index.second ];
          };
          initBlocks( blockAccess, domainEntity, rangeEntity );
        }
      }

      template <class BlockAccess>
      void initBlocks( BlockAccess& blockAccess, const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        auto functor = [ this, &blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          matrices_[ local.first ][ local.second ] = &blockAccess( index );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }

    private:
      // check whether given (row,col) pair is valid
      void check(int localRow, int localCol) const
      {
#ifndef NDEBUG
        const std::size_t row = (int) localRow / littleRows;
        const std::size_t col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        assert( row < matrices_.size() ) ;
        assert( col < matrices_[row].size() );
        assert( lRow < littleRows );
        assert( lCol < littleCols );
#endif
      }

      DofType& getValue(const int localRow, const int localCol)
      {
        const int row = (int) localRow / littleRows;
        const int col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        return (*matrices_[row][col])[lRow][lCol];
      }

    public:
      const DofType get(const int localRow, const int localCol) const
      {
        const int row = (int) localRow / littleRows;
        const int col = (int) localCol / littleCols;
        const int lRow = localRow%littleRows;
        const int lCol = localCol%littleCols;
        return (*matrices_[row][col])[lRow][lCol];
      }

      void scale (const DofType& scalar)
      {
        const size_type nrows = matrices_.size();
        for(size_type i=0; i<nrows; ++i)
        {
          const size_type ncols = matrices_[i].size();
          for(size_type j=0; j<ncols; ++j)
          {
            (*matrices_[i][j]) *= scalar;
          }
        }
      }

      void add(const int localRow, const int localCol , const DofType value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        getValue(localRow,localCol) += value;
      }

      void set(const int localRow, const int localCol , const DofType value)
      {
#ifndef NDEBUG
        check(localRow,localCol);
#endif
        getValue(localRow,localCol) = value;
      }

      //! make unit row (all zero, diagonal entry 1.0 )
      void unitRow(const int localRow)
      {
        const int row = (int) localRow / littleRows;
        const int lRow = localRow%littleRows;

        // clear row
        doClearRow( row, lRow );

        // set diagonal entry to 1
        (*matrices_[row][row])[lRow][lRow] = 1;
      }

      //! clear all entries belonging to local matrix
      void clear ()
      {
        const size_type nrows = matrices_.size();
        for(size_type i=0; i<nrows; ++i)
        {
          const size_type ncols = matrices_[i].size();
          for(size_type j=0; j<ncols; ++j)
          {
            (*matrices_[i][j]) = (DofType) 0;
          }
        }
      }

      //! set matrix row to zero
      void clearRow ( const int localRow )
      {
        const int row = (int) localRow / littleRows;
        const int lRow = localRow%littleRows;

        // clear the row
        doClearRow( row, lRow );
      }

      //! empty as the little matrices are already sorted
      void resort ()
      {}

    protected:
      //! set matrix row to zero
      void doClearRow ( const int row, const int lRow )
      {
        // get number of columns
        const auto col = this->columns();
        for(auto localCol=0; localCol<col; ++localCol)
        {
          const int col = (int) localCol / littleCols;
          const int lCol = localCol%littleCols;
          (*matrices_[row][col])[lRow][lCol] = 0;
        }
      }

    private:
      // special mapper omitting block size
      const RowMapperType& rowMapper_;
      const ColMapperType& colMapper_;

      // number of local matrices
      int numRows_;
      int numCols_;

      // dynamic matrix with pointers to block matrices
      std::vector< std::vector< LittleBlockType* > > matrices_;

      // matrix to build
      const MatrixObjectType& matrixObj_;
    };


    //! MatrixObject handling an istl matrix
    template <class DomainSpaceImp, class RangeSpaceImp, class DomainBlock, class RangeBlock >
    class ISTLMatrixObject
    {
    public:
      //! type of space defining row structure
      typedef DomainSpaceImp DomainSpaceType;
      //! type of space defining column structure
      typedef RangeSpaceImp RangeSpaceType;

      //! type of this pointer
      typedef ISTLMatrixObject< DomainSpaceImp, RangeSpaceImp, DomainBlock, RangeBlock > ThisType;
      typedef ThisType  PreconditionMatrixType;

      typedef typename DomainSpaceType::GridType GridType;

      typedef typename RangeSpaceType  :: EntityType  RangeEntityType ;
      typedef typename DomainSpaceType :: EntityType  DomainEntityType ;

      enum { littleCols = DomainSpaceType :: localBlockSize };
      enum { littleRows = RangeSpaceType  :: localBlockSize };

      typedef FieldMatrix<typename DomainSpaceType :: RangeFieldType, littleRows, littleCols> LittleBlockType;
      typedef LittleBlockType  block_type;
      typedef LittleBlockType  MatrixBlockType;

      typedef ISTLBlockVectorDiscreteFunction< RangeSpaceType, RangeBlock >     RowDiscreteFunctionType;
      typedef ISTLBlockVectorDiscreteFunction< DomainSpaceType, DomainBlock >   ColumnDiscreteFunctionType;

      typedef typename RowDiscreteFunctionType :: DofStorageType    RowBlockVectorType;
      typedef typename ColumnDiscreteFunctionType :: DofStorageType ColumnBlockVectorType;

      typedef typename RangeSpaceType :: BlockMapperType  RowMapperType;
      typedef typename DomainSpaceType :: BlockMapperType ColMapperType;

      //! type of used matrix
      typedef ImprovedBCRSMatrix< LittleBlockType , ColumnDiscreteFunctionType , RowDiscreteFunctionType > MatrixType;
      typedef ISTLParallelMatrixAdapterInterface< MatrixType >    MatrixAdapterType;

      //! type of local matrix
      typedef ISTLLocalMatrix<ThisType> ObjectType;
      friend class ISTLLocalMatrix<ThisType>;

      typedef ThisType LocalMatrixFactoryType;
      typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
      //! type of local matrix
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;
      typedef ColumnObject< ThisType > LocalColumnObjectType;

    protected:
      const DomainSpaceType& domainSpace_;
      const RangeSpaceType&  rangeSpace_;

      // special row mapper
      RowMapperType& rowMapper_;
      // special col mapper
      ColMapperType& colMapper_;

      int size_;
      int sequence_;

      mutable std::unique_ptr< MatrixType > matrix_;

      mutable LocalMatrixStackType localMatrixStack_;

      mutable std::unique_ptr< MatrixAdapterType >     matrixAdap_;
      // overflow fraction for implicit build mode
      const double overflowFraction_;
      const bool threading_ ;

      ISTLSolverParameter param_;
    public:
      ISTLMatrixObject(const ISTLMatrixObject&) = delete;

      //! constructor
      //! \param rowSpace space defining row structure
      //! \param colSpace space defining column structure
      ISTLMatrixObject ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace, const ISTLSolverParameter& param = ISTLSolverParameter() ) :
        domainSpace_(domainSpace)
        , rangeSpace_(rangeSpace)
        , rowMapper_( rangeSpace.blockMapper() )
        , colMapper_( domainSpace.blockMapper() )
        , size_(-1)
        , sequence_(-1)
        , localMatrixStack_( *this )
        , overflowFraction_( param.overflowFraction() )
        , threading_( param.threading() )
        , param_( param )
      {}

    protected:
      //! return reference to system matrix (may not be finalized)
      //! For finalized matrix use exportMatrix
      MatrixType& matrix() const
      {
        assert( matrix_ );
        return *matrix_;
      }

    public:
      bool threading() const { return threading_; }

      void printTexInfo(std::ostream& out) const
      {
        out << "ISTL MatrixObj: ";
        out  << "\\\\ \n";
      }

      //! return matrix adapter object
      std::string preconditionName() const
      {
        return "";
      }

      void createMatrixAdapter ( const ISTLSolverParameter& param ) const
      {
        if( !matrixAdap_ )
        {
          matrixAdap_ = ISTLMatrixAdapterFactory< ThisType > :: matrixAdapter( *this, param );
        }
        assert( matrixAdap_ );
      }

      //! return matrix adapter object
      MatrixAdapterType& matrixAdapter( const ISTLSolverParameter& parameter ) const
      {
        const ISTLSolverParameter* param = dynamic_cast< const ISTLSolverParameter* > (&parameter);
        if( param )
          createMatrixAdapter( *param );
        else
        {
          ISTLSolverParameter newparam( parameter );
          createMatrixAdapter( newparam );
        }

        finalizeAssembly();
        return *matrixAdap_;
      }

      //! return matrix adapter object
      MatrixAdapterType& matrixAdapter() const
      {
        return matrixAdapter( param_ );
      }

    public:
      MatrixType& exportMatrix() const
      {
        finalizeAssembly();
        return matrix();
      }

      bool implicitModeActive() const
      {
        // implicit build mode is only active when the
        // build mode of the matrix is implicit and the
        // matrix is currently being build

        if( matrix().buildMode() == MatrixType::implicit && matrix().buildStage() == MatrixType::building )
          return true;
        else
          return false;
      }

      void finalizeAssembly() const
      {
        // finalize build mode
        const_cast< ThisType& > (*this).compress();
      }

      // compress matrix if not already done before and only in implicit build mode
      void compress( )
      {
        if( implicitModeActive() )
        {
          matrix().compress();
        }
      }

      //! set all matrix entries to zero
      void clear()
      {
        matrix().clear();
        // clean matrix adapter and other helper classes
        removeObj();
      }

      void unitRow( const size_t row )
      {
        matrix().unitRow( row );
      }

    protected:
      template <class Container> // could be std::set or std::vector
      void setUnitRowImpl( const Container &rows, const double diagonal )
      {
        for (const auto& r : rows)
        {
          const std::size_t blockRow( r/(LittleBlockType :: rows) );
          const std::size_t localRowIdx( r%(LittleBlockType :: rows) );
          auto& row = matrix()[blockRow];
          const auto endcol = row.end();
#ifndef NDEBUG
          bool set = false;
#endif
          for (auto col=row.begin(); col!=endcol; ++col)
          {
            for (auto& entry : (*col)[localRowIdx])
              entry = 0;
            if (col.index() == blockRow)
            {
              (*col)[localRowIdx][localRowIdx] = diagonal;
#ifndef NDEBUG
              set = true;
#endif
            }
          }
          assert(set);
        }
      }

    public:
      template <class Container>
      void setUnitRows( const Container& unitRows, const Container& auxRows )
      {
        setUnitRowImpl( unitRows, 1.0 );
        setUnitRowImpl( auxRows,  0.0 );
      }

      //! default reserve method setting implicit build mode
      void reserve()
      {
        reserve( Stencil<DomainSpaceType,RangeSpaceType>(domainSpace_, rangeSpace_), true );
      }

      //! reserve memory for assemble based on the provided stencil
      template <class Set>
      void reserve (const std::vector< Set >& sparsityPattern )
      {
        reserve( StencilWrapper< DomainSpaceType,RangeSpaceType, Set >( sparsityPattern ), false );
      }

      //! reserve memory for assemble based on the provided stencil
      template <class Stencil>
      void reserve(const Stencil &stencil, const bool proposeImplicit = true )
      {
        // if grid sequence number changed, rebuild matrix
        if(sequence_ != domainSpace().sequence())
        {
          removeObj();

          // do not use implicit build mode when multi threading is enabled
          const bool implicit = proposeImplicit && MPIManager::numThreads() == 1;
          if( implicit )
          {
            auto nnz = stencil.maxNonZerosEstimate();
            if( nnz == 0 )
            {
              Stencil tmpStencil( stencil );
              tmpStencil.fill( *(domainSpace_.begin()), *(rangeSpace_.begin()) );
              nnz = tmpStencil.maxNonZerosEstimate();
            }
            matrix_.reset( new MatrixType( rowMapper_.size(), colMapper_.size(), nnz, overflowFraction_ ) );
          }
          else
          {
            matrix_.reset( new MatrixType( rowMapper_.size(), colMapper_.size() ) );
            matrix().createEntries( stencil.globalStencil() );
          }

          sequence_ = domainSpace().sequence();
        }
      }

      //! setup new matrix with hanging nodes taken into account
      template <class HangingNodesType>
      void changeHangingNodes(const HangingNodesType& hangingNodes)
      {
        // create new matrix
        MatrixType* newMatrix = new MatrixType(rowMapper_.size(), colMapper_.size());

        // setup with hanging rows
        newMatrix->setup( *matrix_ , hangingNodes );

        // remove old matrix
        removeObj();

        // store new matrix
        matrix_.reset( newMatrix );
      }

      //! extract diagonal entries of the matrix to a discrete function of type
      //! BlockVectorDiscreteFunction
      void extractDiagonal( ColumnDiscreteFunctionType& diag ) const
      {
        // extract diagonal entries
        matrix().extractDiagonal( diag.blockVector() );
      }

      //! we only have right precondition
      bool rightPrecondition() const
      {
        return false;
      }

      //! apply with discrete functions
      void apply(const ColumnDiscreteFunctionType& arg, RowDiscreteFunctionType& dest) const
      {
        matrixAdapter().apply( arg.blockVector(), dest.blockVector() );
      }

      //! apply with arbitrary discrete functions
      template <class CDF, class RDF>
      void apply(const CDF& arg, RDF& dest) const
      {
        // this can happen when different storages are used for discrete
        // function and matrix/solver objects
        DUNE_THROW(NotImplemented,"ISTLMatrixObj::apply called for DiscreteFunctions not specified in the template list");
      }

      //! resort row numbering in matrix to have ascending numbering
      void resort()
      {}

      //! create precondition matrix does nothing because preconditioner is
      //! created only when requested
      void createPreconditionMatrix()
      {}

      //! print matrix
      void print(std::ostream & s) const
      {
        matrix().print(s);
      }

      const DomainSpaceType& domainSpace() const
      {
        return domainSpace_;
      }
      const RangeSpaceType&  rangeSpace() const
      {
        return rangeSpace_;
      }

      const RowMapperType& rowMapper() const
      {
        return rowMapper_;
      }
      const ColMapperType& colMapper() const
      {
        return colMapper_;
      }

      //! interface method from LocalMatrixFactory
      ObjectType* newObject() const
      {
        return new ObjectType(*this, domainSpace(), rangeSpace());
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

      LocalColumnObjectType localColumn( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType ( *this, domainEntity );
      }

    protected:
      template< class LocalBlock, class Operation >
      void applyToBlock ( const size_t row, const size_t col,
                          const LocalBlock &localBlock,
                          Operation& operation )
      {
        LittleBlockType& block = ( implicitModeActive() ) ? matrix().entry( row, col ) : matrix()[ row ][ col ];
        for( int i  = 0; i < littleRows; ++i )
          for( int j = 0; j < littleCols; ++j )
            operation( block[ i ][ j ], localBlock[ i ][ j ] );
      }

      template< class LocalBlock >
      void setBlock ( const size_t row, const size_t col,
                      const LocalBlock &localBlock )
      {
        typedef typename DomainSpaceType :: RangeFieldType Field;
        auto copy = [] ( Field& a, const typename LocalBlock::field_type& b ) { a = b; };
        applyToBlock( row, col, localBlock, copy );
      }

      template< class LocalBlock >
      void addBlock ( const size_t row, const size_t col,
                      const LocalBlock &localBlock )
      {
        typedef typename DomainSpaceType :: RangeFieldType Field;
        auto add = [] ( Field& a, const typename LocalBlock::field_type& b ) { a += b; };
        applyToBlock( row, col, localBlock, add );
      }

    public:
      template< class LocalMatrix, class Operation >
      void applyToLocalMatrix ( const DomainEntityType &domainEntity,
                                const RangeEntityType &rangeEntity,
                                const LocalMatrix &localMat,
                                Operation& operation )
      {
        typedef typename MatrixType::size_type Index;
        if( implicitModeActive() )
        {
          auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
          {
            return matrix().entry( index.first, index.second );
          };
          auto functor = [ &localMat, blockAccess, operation ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
          {
            LittleBlockType& block = blockAccess( index );
            for( int i  = 0; i < littleRows; ++i )
              for( int j = 0; j < littleCols; ++j )
                operation( block[ i ][ j ], localMat.get( local.first * littleRows + i, local.second *littleCols + j ) );
          };
          rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
        }
        else
        {
          auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
          {
            return matrix()[ index.first][ index.second ];
          };
          auto functor = [ &localMat, blockAccess, operation ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
          {
            LittleBlockType& block = blockAccess( index );
            for( int i  = 0; i < littleRows; ++i )
              for( int j = 0; j < littleCols; ++j )
                operation( block[ i ][ j ], localMat.get( local.first * littleRows + i, local.second *littleCols + j ) );
          };
          rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
        }
      }

      template< class LocalMatrix >
      void setLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        typedef typename DomainSpaceType :: RangeFieldType  RangeFieldType;
        auto operation = [] ( RangeFieldType& a, const RangeFieldType& b )
        {
          a = b;
        };
        applyToLocalMatrix( domainEntity, rangeEntity, localMat, operation );
      }

      template< class LocalMatrix >
      void addLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        typedef typename DomainSpaceType :: RangeFieldType  RangeFieldType;
        auto operation = [] ( RangeFieldType& a, const RangeFieldType& b )
        {
          a += b;
        };
        applyToLocalMatrix( domainEntity, rangeEntity, localMat, operation );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat, const Scalar &s  )
      {
        typedef typename DomainSpaceType :: RangeFieldType  RangeFieldType;
        auto operation = [ &s ] ( RangeFieldType& a, const RangeFieldType& b )
        {
          a += s * b;
        };
        applyToLocalMatrix( domainEntity, rangeEntity, localMat, operation );
      }

      template< class LocalMatrix >
      void getLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, LocalMatrix &localMat ) const
      {
        // make sure that matrix is in compressed state if build mode is implicit
        finalizeAssembly();

        typedef typename MatrixType::size_type Index;
        auto functor = [ &localMat, this ] ( std::pair< int, int > local, const std::pair< Index, Index > &global )
        {
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              localMat.set( local.first * littleRows + i, local.second *littleCols + j,
                matrix()[ global.first ][ global.second ][i][j] );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }

    protected:
      void preConErrorMsg(int preCon) const
      {
        exit(1);
      }

      void removeObj ()
      {
        matrixAdap_.reset( nullptr );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLMATRIXWRAPPER_HH
