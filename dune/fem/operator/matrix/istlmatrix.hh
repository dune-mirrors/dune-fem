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

//- Dune istl includes
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/preconditioners.hh>

//- Dune fem includes
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/preconditionerwrapper.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/storage/objectstack.hh>

#include <dune/fem/operator/matrix/istlmatrixadapter.hh>
#include <dune/fem/operator/matrix/functor.hh>

namespace Dune
{

  namespace Fem
  {


    //! forward declrations
    template< class MatrixObject >
    class ISTLLocalMatrix;

    template< class RowSpaceImp, class ColSpaceImp >
    class ISTLMatrixObject;

    struct ISTLMatrixParameter
      : public MatrixParameter
    {

      ISTLMatrixParameter( const std::string keyPrefix = "istl." )
        : keyPrefix_( keyPrefix )
      {}

      virtual double overflowFraction () const
      {
        return Parameter::getValue< double >( keyPrefix_ + "matrix.overflowfraction", 1.0 );
      }

      virtual int numIterations () const
      {
        return Parameter::getValue< int >( keyPrefix_ + "preconditioning.iterations", 5 );
      }

      virtual double relaxation () const
      {
        return Parameter::getValue< double >( keyPrefix_ + "preconditioning.relaxation", 1.1 );
      }

      virtual int method () const
      {
        static const std::string preConTable[]
          = { "none", "ssor", "sor", "ilu-0", "ilu-n", "gauss-seidel", "jacobi", "amg-ilu-0", "amg-ilu-n", "amg-jacobi" };
        return Parameter::getEnum(  keyPrefix_ + "preconditioning.method", preConTable, 0 );
      }

      virtual std::string preconditionName() const
      {
        static const std::string preConTable[]
          = { "None", "SSOR", "SOR", "ILU-0", "ILU-n", "Gauss-Seidel", "Jacobi", "AMG-ILU-0", "AMG-ILU-n", "AMG-Jacobi" };
        const int precond = method();
        std::stringstream tmp;
        tmp << preConTable[precond];

        if( precond != 3 )
          tmp << " n=" << numIterations();
        tmp << " relax=" << relaxation();
        return tmp.str();
      }

     private:
      std::string keyPrefix_;

    };

    ///////////////////////////////////////////////////////
    // --ISTLMatrixHandle
    //////////////////////////////////////////////////////
    template <class LittleBlockType, class RowDiscreteFunctionImp, class ColDiscreteFunctionImp = RowDiscreteFunctionImp>
    class ImprovedBCRSMatrix : public BCRSMatrix<LittleBlockType>
    {
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

        //! increment block level counter
        enum {
          //! The number of blocklevels the matrix contains.
          blocklevel = BaseType :: blocklevel
        };

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
        typedef typename RangeSpaceType :: GridType :: Traits :: CollectiveCommunication   CollectiveCommunictionType ;

        typedef typename BaseType :: BuildMode BuildMode ;

      public:
        //! constructor used by ISTLMatrixObject to build matrix in implicit mode
        ImprovedBCRSMatrix(size_type rows, size_type cols, size_type nnz, double overflowFraction) :
          BaseType (rows, cols, nnz, overflowFraction, BaseType::implicit)
        {}

        //! constuctor using old row_wise assembly (and used by ILU preconditioner)
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

        template <class RowKeyType, class ColKeyType>
        void createEntries(const std::map<RowKeyType , std::set<ColKeyType> >& indices)
        {
          // not insert map of indices into matrix
          auto endcreate = this->createend();
          for(auto create = this->createbegin(); create != endcreate; ++create)
          {
            const auto it = indices.find( create.index() );
            if (it == indices.end() )
              continue;
            const auto& localIndices = it->second;
            const auto end = localIndices.end();
            // insert all indices for this row
            for (auto it = localIndices.begin(); it != end; ++it)
              create.insert( *it );
          }
        }

        //! clear Matrix, i.e. set all entires to 0
        void clear()
        {
          for (auto& row : *this)
            for (auto& entry : row)
              entry = 0;
        }

        //! setup like the old matrix but remove rows with hanging nodes
        template <class HangingNodesType>
        void setup(ThisType& oldMatrix, const HangingNodesType& hangingNodes)
        {
          // necessary because element traversal not necessaryly is in ascending order
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

      ISTLLocalMatrix ( const MatrixObjectType& mObj, const DomainSpaceType& domainSpace, const RangeSpaceType& rangeSpace )
        : BaseType( domainSpace, rangeSpace ),
          rowMapper_( rangeSpace.blockMapper() ),
          colMapper_( domainSpace.blockMapper() ),
          numRows_( rowMapper_.maxNumDofs() ),
          numCols_( colMapper_.maxNumDofs() ),
          matrixObj_( mObj )
      {}

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
        matrices_.resize( numRows_, numCols_, nullptr );

        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( matrixObj_.implicitModeActive() )
            return matrixObj_.matrix().entry( index.first, index.second );
          else
            return matrixObj_.matrix()[ index.first ][  index.second ];
        };

        auto functor = [ this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
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
        for(auto i=0; i<matrices_.size(); ++i)
          for(auto j=0; j<matrices_[i].size(); ++j)
            (*matrices_[i][j]) *= scalar;
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
        for(auto i=0; i<matrices_.size(); ++i)
          for(auto j=0; j<matrices_[i].size(); ++j)
            (*matrices_[i][j]) = (DofType) 0;
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
      // special mapper omiting block size
      const RowMapperType& rowMapper_;
      const ColMapperType& colMapper_;

      // number of local matrices
      int numRows_;
      int numCols_;

      // dynamic matrix with pointers to block matrices
      Dune::DynamicMatrix< LittleBlockType* > matrices_;

      // matrix to build
      const MatrixObjectType& matrixObj_;
    };



    template <class RowSpaceImp, class ColSpaceImp = RowSpaceImp>
    struct ISTLMatrixTraits
    {
      typedef RowSpaceImp RangeSpaceType;
      typedef ColSpaceImp DomainSpaceType;
      typedef ISTLMatrixTraits<DomainSpaceType,RangeSpaceType> ThisType;

      typedef ISTLMatrixObject<DomainSpaceType,RangeSpaceType> MatrixObjectType;
    };

    //! MatrixObject handling an istl matrix
    template <class DomainSpaceImp, class RangeSpaceImp>
    class ISTLMatrixObject
    {
    public:
      //! type of space defining row structure
      typedef DomainSpaceImp DomainSpaceType;
      //! type of space defining column structure
      typedef RangeSpaceImp RangeSpaceType;

      //! type of this pointer
      typedef ISTLMatrixObject<DomainSpaceType,RangeSpaceType> ThisType;

      typedef typename DomainSpaceType::GridType GridType;

      typedef typename RangeSpaceType :: EntityType  RangeEntityType ;
      typedef typename DomainSpaceType :: EntityType DomainEntityType ;

      enum { littleCols = DomainSpaceType :: localBlockSize };
      enum { littleRows = RangeSpaceType :: localBlockSize };

      typedef FieldMatrix<typename DomainSpaceType :: RangeFieldType, littleRows, littleCols> LittleBlockType;

      typedef ISTLBlockVectorDiscreteFunction< RangeSpaceType >      RowDiscreteFunctionType;
      typedef ISTLBlockVectorDiscreteFunction< DomainSpaceType >     ColumnDiscreteFunctionType;

    protected:
      typedef typename RowDiscreteFunctionType :: DofStorageType    RowBlockVectorType;
      typedef typename ColumnDiscreteFunctionType :: DofStorageType ColumnBlockVectorType;

      typedef typename RangeSpaceType :: BlockMapperType  RowMapperType;
      typedef typename DomainSpaceType :: BlockMapperType ColMapperType;

    public:
      //! type of used matrix
      typedef ImprovedBCRSMatrix< LittleBlockType , ColumnDiscreteFunctionType , RowDiscreteFunctionType > MatrixType;
      typedef typename ISTLParallelMatrixAdapter<MatrixType,RangeSpaceType>::Type MatrixAdapterType;
      typedef typename MatrixAdapterType :: ParallelScalarProductType ParallelScalarProductType;

      //! type of local matrix
      typedef ISTLLocalMatrix<ThisType> ObjectType;
      typedef ThisType LocalMatrixFactoryType;
      typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
      //! type of local matrix
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;
      typedef ColumnObject< ThisType > LocalColumnObjectType;

    protected:
      const DomainSpaceType & domainSpace_;
      const RangeSpaceType & rangeSpace_;

      // sepcial row mapper
      RowMapperType& rowMapper_;
      // special col mapper
      ColMapperType& colMapper_;

      int size_;

      int sequence_;

      mutable std::unique_ptr< MatrixType > matrix_;

      mutable LocalMatrixStackType localMatrixStack_;

      mutable std::unique_ptr< MatrixAdapterType > matrixAdap_;
      mutable std::unique_ptr< ColumnBlockVectorType > Arg_;
      mutable std::unique_ptr< RowBlockVectorType >    Dest_;
      // overflow fraction for implicit build mode
      const double overflowFraction_;

    public:
      ISTLMatrixObject(const ISTLMatrixObject&) = delete;

      //! constructor
      //! \param rowSpace space defining row structure
      //! \param colSpace space defining column structure
      ISTLMatrixObject ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace, const MatrixParameter& param = ISTLMatrixParameter() ) :
        domainSpace_(domainSpace)
        , rangeSpace_(rangeSpace)
        , rowMapper_( rangeSpace.blockMapper() )
        , colMapper_( domainSpace.blockMapper() )
        , size_(-1)
        , sequence_(-1)
        , localMatrixStack_( *this )
        , overflowFraction_( param.overflowFraction() )
      {}

      ThisType &systemMatrix () { return *this; }
      const ThisType &systemMatrix () const { return *this; }

      //! return reference to system matrix
      MatrixType & matrix() const
      {
        assert( matrix_ );
        return *matrix_;
      }

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

      void createMatrixAdapter () const
      {
        if( !matrixAdap_ )
          matrixAdap_.reset( new MatrixAdapterType( matrixAdapterObject() ) );
        assert( matrixAdap_ );
      }

      //! return matrix adapter object
      const MatrixAdapterType& matrixAdapter() const
      {
        createMatrixAdapter();
        return *matrixAdap_;
      }
      MatrixAdapterType& matrixAdapter()
      {
        createMatrixAdapter();
        return *matrixAdap_;
      }

    protected:
      MatrixAdapterType matrixAdapterObject() const
      {
        // need some precondition object - empty here
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        return MatrixAdapterType( matrix(), domainSpace(), rangeSpace(), PreConType() );
      }

    public:
      bool implicitModeActive() const
      {
        // implicit build mode is only active when the
        // build mode of the matrix is implicit and the
        // matrix is currently being build

        if( matrix().buildMode()  == MatrixType::implicit && matrix().buildStage() == MatrixType::building )
          return true;
        else
          return false;
      }

      // compress matrix if not already done before and only in implicit build mode
      void communicate( )
      {
        if( implicitModeActive() )
          matrix().compress();
      }

      //! return true, because in case of no preconditioning we have empty
      //! preconditioner (used by OEM methods)
      bool hasPreconditionMatrix() const
      {
        return true;
      }

      //! set all matrix entries to zero
      void clear()
      {
        matrix().clear();
        // clean matrix adapter and other helper classes
        removeObj();
      }

      //! reserve memory for assemble based on the provided stencil
      template <class Stencil>
      void reserve(const Stencil &stencil, const bool implicit = true )
      {
        // if grid sequence number changed, rebuild matrix
        if(sequence_ != domainSpace().sequence())
        {
          removeObj();

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

      //! mult method for OEM Solver
      void multOEM(const double* arg, double* dest) const
      {
        createBlockVectors();

        RowBlockVectorType& Arg = *Arg_;
        ColumnBlockVectorType &Dest = *Dest_;

        std::copy_n( arg, Arg.size(), Arg.dbegin() );

        // call mult of matrix adapter
        matrixAdapter().apply( Arg, Dest );

        std::copy( Dest.dbegin(), Dest.dend(), dest );
      }

      //! apply with discrete functions
      void apply(const ColumnDiscreteFunctionType& arg, RowDiscreteFunctionType& dest) const
      {
        matrixAdapter().apply( arg.blockVector(), dest.blockVector() );
      }

      //! apply with arbitrary discrete functions calls multOEM
      template <class RowDFType, class ColDFType>
      void apply(const RowDFType& arg, ColDFType& dest) const
      {
        multOEM( arg.leakPointer(), dest.leakPointer ());
      }

      //! mult method of matrix object used by oem solver
      template <class ColumnLeakPointerType, class RowLeakPointerType>
      void multOEM(const ColumnLeakPointerType& arg, RowLeakPointerType& dest) const
      {
        DUNE_THROW(NotImplemented,"Method has been removed");
      }

      //! dot method for OEM Solver
      double ddotOEM(const double* v, const double* w) const
      {
        createBlockVectors();

        RowBlockVectorType&    V = *Arg_;
        ColumnBlockVectorType& W = *Dest_;

        std::copy_n( v, V.size(), V.dbegin() );
        std::copy_n( w, W.size(), W.dbegin() );

#if HAVE_MPI
        // in parallel use scalar product of discrete functions
        ISTLBlockVectorDiscreteFunction< DomainSpaceType > vF("ddotOEM:vF", domainSpace(), V );
        ISTLBlockVectorDiscreteFunction< RangeSpaceType  > wF("ddotOEM:wF", rangeSpace(), W );
        return vF.scalarProductDofs( wF );
#else
        return V * W;
#endif
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

      //! return local matrix object
      LocalMatrixType localMatrix( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity ) const
      {
        return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
      }
      LocalColumnObjectType localColumn( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType ( *this, domainEntity );
      }

      template< class LocalMatrix >
      void addLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( implicitModeActive() )
            return matrix().entry( index.first, index.second );
          else
            return matrix()[ index.first ][  index.second ];
        };

        auto functor = [ &localMat, this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          LittleBlockType& block = blockAccess( index );
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              block[ i ][ j ] += localMat.get( local.first * littleRows + i, local.second *littleCols + j );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }


      template< class LocalMatrix >
      void setLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( implicitModeActive() )
            return matrix().entry( index.first, index.second );
          else
            return matrix()[ index.first ][  index.second ];
        };

        auto functor = [ &localMat, this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          LittleBlockType& block = blockAccess( index );
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              block[ i ][ j ] = localMat.get( local.first * littleRows + i, local.second *littleCols + j );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat, const Scalar &s  )
      {
        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( implicitModeActive() )
            return matrix().entry( index.first, index.second );
          else
            return matrix()[ index.first ][  index.second ];
        };

        auto functor = [ &localMat, &s, this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          LittleBlockType& block = blockAccess( index );
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              block[ i ][ j ] += s * localMat.get( local.first * littleRows + i, local.second *littleCols + j );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix >
      void getLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, LocalMatrix &localMat ) const
      {
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
        Dest_.reset( nullptr );
        Arg_.reset( nullptr );
        matrixAdap_.reset( nullptr );
      }

      void createBlockVectors () const
      {
        if( !Arg_  )
          Arg_.reset( new RowBlockVectorType( rowMapper_.size() ) );
        if( !Dest_ )
          Dest_.reset( new ColumnBlockVectorType( colMapper_.size() ) );

        createMatrixAdapter ();
      }
    };



    //! MatrixObject handling an istl matrix
    template <class SpaceImp>
    class ISTLMatrixObject<SpaceImp,SpaceImp>
    {
    public:
      typedef SpaceImp DomainSpaceImp;
      typedef SpaceImp RangeSpaceImp;
      //! type of space defining row structure
      typedef DomainSpaceImp DomainSpaceType;
      //! type of space defining column structure
      typedef RangeSpaceImp RangeSpaceType;

      //! type of this pointer
      typedef ISTLMatrixObject<DomainSpaceType,RangeSpaceType> ThisType;

      typedef typename DomainSpaceType::GridType GridType;

      typedef typename RangeSpaceType :: EntityType RangeEntityType ;
      typedef typename DomainSpaceType :: EntityType DomainEntityType;

      enum { littleRows = DomainSpaceType :: localBlockSize };
      enum { littleCols = RangeSpaceType :: localBlockSize };

      typedef FieldMatrix<typename DomainSpaceType :: RangeFieldType, littleRows, littleCols> LittleBlockType;

      typedef ISTLBlockVectorDiscreteFunction< DomainSpaceType >     RowDiscreteFunctionType;
      typedef ISTLBlockVectorDiscreteFunction< RangeSpaceType >  ColumnDiscreteFunctionType;

    protected:
      typedef typename RowDiscreteFunctionType    :: DofContainerType  RowBlockVectorType;
      typedef typename ColumnDiscreteFunctionType :: DofContainerType  ColumnBlockVectorType;

      typedef typename DomainSpaceType :: BlockMapperType RowMapperType;
      typedef typename RangeSpaceType :: BlockMapperType ColMapperType;

    public:
      //! type of used matrix
      typedef ImprovedBCRSMatrix< LittleBlockType , RowDiscreteFunctionType , ColumnDiscreteFunctionType > MatrixType;
      typedef typename ISTLParallelMatrixAdapter<MatrixType,DomainSpaceType>::Type MatrixAdapterType;
      // get preconditioner type from MatrixAdapterType
      typedef ThisType PreconditionMatrixType;
      typedef typename MatrixAdapterType :: ParallelScalarProductType ParallelScalarProductType;

      //! type of local matrix
      typedef ISTLLocalMatrix<ThisType> ObjectType;
      typedef ThisType LocalMatrixFactoryType;
      typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;
      //! type of local matrix
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;
      typedef ColumnObject< ThisType > LocalColumnObjectType;

    protected:
      const DomainSpaceType & domainSpace_;
      const RangeSpaceType & rangeSpace_;

      // sepcial row mapper
      RowMapperType& rowMapper_;
      // special col mapper
      ColMapperType& colMapper_;

      int size_;

      int sequence_;

      mutable std::unique_ptr< MatrixType > matrix_;

      ParallelScalarProductType scp_;

      int numIterations_;
      double relaxFactor_;

      enum ISTLPreConder_Id { none  = 0 ,      // no preconditioner
                              ssor  = 1 ,      // SSOR preconditioner
                              sor   = 2 ,      // SOR preconditioner
                              ilu_0 = 3 ,      // ILU-0 preconditioner
                              ilu_n = 4 ,      // ILU-n preconditioner
                              gauss_seidel= 5, // Gauss-Seidel preconditioner
                              jacobi = 6,      // Jacobi preconditioner
                              amg_ilu_0 = 7,   // AMG with ILU-0 smoother
                              amg_ilu_n = 8,   // AMG with ILU-n smoother
                              amg_jacobi = 9   // AMG with Jacobi smoother
      };

      ISTLPreConder_Id preconditioning_;

      mutable LocalMatrixStackType localMatrixStack_;

      mutable std::unique_ptr< MatrixAdapterType > matrixAdap_;
      mutable std::unique_ptr< RowBlockVectorType > Arg_;
      mutable std::unique_ptr< ColumnBlockVectorType > Dest_;
      // overflow fraction for implicit build mode
      const double overflowFraction_;
      const MatrixParameter& param_;

    public:
      ISTLMatrixObject(const ISTLMatrixObject&) = delete;

      //! constructor
      //! \param rowSpace space defining row structure
      //! \param colSpace space defining column structure
      //! \param param istl matrix parameters for preconditioning
      //!         - Preconditioning: {0,1,2,3,4,5,6} put -1 to get info
      //!         - Pre-iteration: number of iteration of preconditioner
      //!         - Pre-relaxation: relaxation factor
      ISTLMatrixObject ( const DomainSpaceType &rowSpace, const RangeSpaceType &colSpace, const MatrixParameter& param = ISTLMatrixParameter() )
        : domainSpace_(rowSpace)
        , rangeSpace_(colSpace)
        // create scp to have at least one instance
        // otherwise instance will be deleted during setup
        // get new mappers with number of dofs without considerung block size
        , rowMapper_( rowSpace.blockMapper() )
        , colMapper_( colSpace.blockMapper() )
        , size_(-1)
        , sequence_(-1)
        , scp_(rangeSpace())
        , numIterations_( param.numIterations() )
        , relaxFactor_( param.relaxation() )
        , preconditioning_( (ISTLPreConder_Id)param.method() )
        , localMatrixStack_( *this )
        , overflowFraction_( param.overflowFraction() )
        , param_( param )
      {
        assert( rowMapper_.size() == colMapper_.size() );
      }

      ThisType &systemMatrix () { return *this; }
      const ThisType &systemMatrix () const { return *this; }

      //! return reference to system matrix
      MatrixType & matrix() const
      {
        assert( matrix_ );
        return *matrix_;
      }

      void printTexInfo(std::ostream& out) const
      {
        out << "ISTL MatrixObj: ";
        out << " preconditioner = " << preconditionName() ;
        out  << "\\\\ \n";
      }

      //! return matrix adapter object
      std::string preconditionName() const
      {
        return param_.preconditionName();
      }

      template <class PreconditionerType>
      MatrixAdapterType createMatrixAdapter(const PreconditionerType* preconditioning, std::size_t numIterations) const
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        PreConType preconAdapter(matrix(), numIterations, relaxFactor_, preconditioning );
        return MatrixAdapterType(matrix(), domainSpace(), rangeSpace(), preconAdapter );
      }

      template <class PreconditionerType>
      MatrixAdapterType createAMGMatrixAdapter(const PreconditionerType* preconditioning, std::size_t numIterations) const
      {
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        PreConType preconAdapter(matrix(), numIterations, relaxFactor_, preconditioning, domainSpace().gridPart().comm() );
        return MatrixAdapterType(matrix(), domainSpace(), rangeSpace(), preconAdapter );
      }

      //! return matrix adapter object
      const MatrixAdapterType& matrixAdapter() const
      {
        if( !matrixAdap_ )
          matrixAdap_.reset( new MatrixAdapterType( matrixAdapterObject() ) );
        return *matrixAdap_;
      }
      MatrixAdapterType& matrixAdapter()
      {
        if( !matrixAdap_ )
          matrixAdap_.reset( new MatrixAdapterType( matrixAdapterObject() ) );
        return *matrixAdap_;
      }

    protected:
      MatrixAdapterType matrixAdapterObject() const
      {
#ifndef DISABLE_ISTL_PRECONDITIONING
        const auto procs = domainSpace().gridPart().comm().size();

        typedef typename MatrixType :: BaseType ISTLMatrixType;
        typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
        // no preconditioner
        if( preconditioning_ == none )
        {
          return MatrixAdapterType(matrix(), domainSpace(),rangeSpace(), PreConType() );
        }
        // SSOR
        else if( preconditioning_ == ssor )
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqSSOR not working in parallel computations");

          typedef SeqSSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // SOR
        else if(preconditioning_ == sor )
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqSOR not working in parallel computations");

          typedef SeqSOR<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // ILU-0
        else if(preconditioning_ == ilu_0)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqILU0 not working in parallel computations");

          typedef FemSeqILU0<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // ILU-n
        else if(preconditioning_ == ilu_n)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqILUn not working in parallel computations");

          typedef SeqILUn<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // Gauss-Seidel
        else if(preconditioning_ == gauss_seidel)
        {
          if( procs > 1 )
            DUNE_THROW(InvalidStateException,"ISTL::SeqGS not working in parallel computations");

          typedef SeqGS<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // Jacobi
        else if(preconditioning_ == jacobi)
        {
          if( numIterations_ == 1 ) // diagonal preconditioning
          {
            typedef FemDiagonalPreconditioner< ThisType, RowBlockVectorType, ColumnBlockVectorType > PreconditionerType;
            typedef typename MatrixAdapterType :: PreconditionAdapterType PreConType;
            PreConType preconAdapter( matrix(), new PreconditionerType( *this ) );
            return MatrixAdapterType( matrix(), domainSpace(), rangeSpace(), preconAdapter );
          }
          else if ( procs == 1 )
          {
            typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
            return createMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
          }
          else
          {
            DUNE_THROW(InvalidStateException,"ISTL::SeqJac(Jacobi) only working with istl.preconditioning.iterations: 1 in parallel computations");
          }
        }
        // AMG ILU-0
        else if(preconditioning_ == amg_ilu_0)
        {
          // use original SeqILU0 because of some AMG traits classes.
          typedef SeqILU0<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createAMGMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // AMG ILU-n
        else if(preconditioning_ == amg_ilu_n)
        {
          typedef SeqILUn<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType> PreconditionerType;
          return createAMGMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        // AMG Jacobi
        else if(preconditioning_ == amg_jacobi)
        {
          typedef SeqJac<ISTLMatrixType,RowBlockVectorType,ColumnBlockVectorType,1> PreconditionerType;
          return createAMGMatrixAdapter( (PreconditionerType*)nullptr, numIterations_ );
        }
        else
        {
          preConErrorMsg(preconditioning_);
        }
#endif

        return MatrixAdapterType(matrix(), domainSpace(), rangeSpace(), PreConType() );
      }

    public:
      bool implicitModeActive() const
      {
        // implicit build mode is only active when the
        // build mode of the matrix is implicit and the
        // matrix is currently being build
        if( matrix().buildMode()  == MatrixType::implicit && matrix().buildStage() == MatrixType::building )
          return true;
        else
          return false;
      }

      // compress matrix if not already done before and only in implicit build mode
      void communicate( )
      {
        if( implicitModeActive() )
          matrix().compress();
      }

      //! return true, because in case of no preconditioning we have empty
      //! preconditioner (used by OEM methods)
      bool hasPreconditionMatrix() const
      {
        return (preconditioning_ != none);
      }

      //! return reference to preconditioner object (used by OEM methods)
      const PreconditionMatrixType& preconditionMatrix() const
      {
        return *this;
      }

      //! set all matrix entries to zero
      void clear()
      {
        matrix().clear();
        // clean matrix adapter and other helper classes
        removeObj();
      }

      //! reserve memory for assemble based on the provided stencil
      template <class Stencil>
      void reserve( const Stencil &stencil, const bool implicit = true )
      {
        // if grid sequence number changed, rebuild matrix
        if(sequence_ != domainSpace().sequence())
        {
          removeObj();

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
        return true;
      }

      //! precondition method for OEM Solvers
      //! not fast but works, double is copied to block vector
      //! and after application copied back
      void precondition(const double* arg, double* dest) const
      {
        createBlockVectors();

        assert( Arg_ );
        assert( Dest_ );

        RowBlockVectorType& Arg = *Arg_;
        ColumnBlockVectorType & Dest = *Dest_;

        std::copy_n( arg, Arg.size(), Arg.dbegin());

        // set Dest to zero
        Dest = 0;

        assert( matrixAdap_ );
        // not parameter swaped for preconditioner
        matrixAdap_->preconditionAdapter().apply(Dest , Arg);

        std::copy( Dest.dbegin(), Dest.dend(), dest );
      }

      //! mult method for OEM Solver
      void multOEM(const double* arg, double* dest) const
      {
        createBlockVectors();

        assert( Arg_ );
        assert( Dest_ );

        RowBlockVectorType& Arg = *Arg_;
        ColumnBlockVectorType & Dest = *Dest_;

        std::copy_n( arg, Arg.size(), Arg.dbegin() );

        // call mult of matrix adapter
        assert( matrixAdap_ );
        matrixAdap_->apply( Arg, Dest );

        std::copy( Dest.dbegin(), Dest.dend(), dest );
      }

      //! apply with discrete functions
      void apply(const RowDiscreteFunctionType& arg, ColumnDiscreteFunctionType& dest) const
      {
        createMatrixAdapter();
        assert( matrixAdap_ );
        matrixAdap_->apply( arg.blockVector(), dest.blockVector() );
      }

      //! apply with arbitrary discrete functions calls multOEM
      template <class RowDFType, class ColDFType>
      void apply(const RowDFType& arg, ColDFType& dest) const
      {
        multOEM( arg.leakPointer(), dest.leakPointer ());
      }

      //! dot method for OEM Solver
      double ddotOEM(const double* v, const double* w) const
      {
        createBlockVectors();

        assert( Arg_ );
        assert( Dest_ );

        RowBlockVectorType&    V = *Arg_;
        ColumnBlockVectorType& W = *Dest_;

        std::copy_n( v, V.size(), V.dbegin() );
        std::copy_n( w, W.size(), W.dbegin() );

#if HAVE_MPI
        // in parallel use scalar product of discrete functions
        ISTLBlockVectorDiscreteFunction< DomainSpaceType > vF("ddotOEM:vF", domainSpace(), V );
        ISTLBlockVectorDiscreteFunction< RangeSpaceType  > wF("ddotOEM:wF", rangeSpace(), W );
        return vF.scalarProductDofs( wF );
#else
        return V * W;
#endif
      }

      //! resort row numbering in matrix to have ascending numbering
      void resort()
      {}

      //! create precondition matrix does nothing because preconditioner is created only when requested
      void createPreconditionMatrix()
      {}

      //! print matrix
      void print(std::ostream & s) const
      {
        matrix().print(s);
      }

      const DomainSpaceType& domainSpace() const { return domainSpace_; }
      const RangeSpaceType&  rangeSpace() const { return rangeSpace_; }

      const RowMapperType& rowMapper() const { return rowMapper_; }
      const ColMapperType& colMapper() const { return colMapper_; }

      //! interface method from LocalMatrixFactory
      ObjectType* newObject() const
      {
        return new ObjectType(*this, domainSpace(), rangeSpace());
      }

      //! return local matrix object
      LocalMatrixType localMatrix(const DomainEntityType& domainEntity, const RangeEntityType &rangeEntity ) const
      {
        return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
      }
      LocalColumnObjectType localColumn( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType ( *this, domainEntity );
      }

      template< class LocalMatrix >
      void addLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( implicitModeActive() )
            return matrix().entry( index.first, index.second );
          else
            return matrix()[ index.first ][  index.second ];
        };

        auto functor = [ &localMat, this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          LittleBlockType& block = blockAccess( index );
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              block[ i ][ j ] += localMat.get( local.first * littleRows + i, local.second *littleCols + j );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }


      template< class LocalMatrix >
      void setLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( implicitModeActive() )
            return matrix().entry( index.first, index.second );
          else
            return matrix()[ index.first ][  index.second ];
        };

        auto functor = [ &localMat, this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          LittleBlockType& block = blockAccess( index );
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              block[ i ][ j ] = localMat.get( local.first * littleRows + i, local.second *littleCols + j );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat, const Scalar &s  )
      {
        typedef typename MatrixType::size_type Index;
        auto blockAccess = [ this ] ( const std::pair< Index, Index > &index ) -> LittleBlockType&
        {
          if( implicitModeActive() )
            return matrix().entry( index.first, index.second );
          else
            return matrix()[ index.first ][  index.second ];
        };

        auto functor = [ &localMat, &s, this, blockAccess ] ( std::pair< int, int > local, const std::pair< Index, Index > &index )
        {
          LittleBlockType& block = blockAccess( index );
          for( std::size_t i  = 0; i < littleRows; ++i )
            for( std::size_t j = 0; j < littleCols; ++j )
              block[ i ][ j ] += s * localMat.get( local.first * littleRows + i, local.second *littleCols + j );
        };

        rowMapper_.mapEach( rangeEntity, makePairFunctor( colMapper_, domainEntity, functor ) );
      }

      template< class LocalMatrix >
      void getLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, LocalMatrix &localMat ) const
      {
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
        std::cerr << "ERROR: Wrong precoditioning number (p = " << preCon;
        std::cerr <<") in ISTLMatrixObject! \n";
        std::cerr <<"Valid values are: \n";
        std::cerr <<"0 == no \n";
        std::cerr <<"1 == SSOR \n";
        std::cerr <<"2 == SOR \n";
        std::cerr <<"3 == ILU-0 \n";
        std::cerr <<"4 == ILU-n \n";
        std::cerr <<"5 == Gauss-Seidel \n";
        std::cerr <<"6 == Jacobi \n";
        assert(false);
        exit(1);
      }

      void removeObj ()
      {
        Dest_.reset( nullptr );
        Arg_.reset( nullptr );
        matrixAdap_.reset( nullptr );
      }

      void createBlockVectors() const
      {
        if( !Arg_ || !Dest_ )
        {
          Arg_.reset( new RowBlockVectorType( rowMapper_.size() ) );
          Dest_.reset( new ColumnBlockVectorType( colMapper_.size() ) );
        }

        createMatrixAdapter ();
      }

      void createMatrixAdapter () const
      {
        if( !matrixAdap_ )
          matrixAdap_.reset( new MatrixAdapterType(matrixAdapterObject()) );
      }

    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_ISTL

#endif // #ifndef DUNE_FEM_ISTLMATRIXWRAPPER_HH
