// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
#define DUNE_FEM_PETSCLINEAROPERATOR_HH

#include <iostream>
#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/fem/misc/functor.hh>
#include <dune/fem/operator/common/localmatrix.hh> 
#include <dune/fem/operator/common/localmatrixwrapper.hh>

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/operator/common/stencil.hh>

#if defined HAVE_PETSC

#include "petscmat.h"


namespace Dune 
{
  namespace Fem 
  {

    /* ========================================
     * class PetscLinearOperator
     */
    template< typename DomainFunction, typename RangeFunction >
    class PetscLinearOperator 
    : public Fem::Operator< DomainFunction, RangeFunction >
    {
      typedef PetscLinearOperator< DomainFunction, RangeFunction > ThisType;
    public:
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;
      typedef typename DomainFunctionType::RangeFieldType DomainFieldType;
      typedef typename RangeFunctionType::RangeFieldType RangeFieldType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

      typedef typename DomainSpaceType::GridPartType::template Codim< 0 >::EntityType RowEntityType;
      typedef typename RangeSpaceType::GridPartType::template Codim< 0 >::EntityType ColumnEntityType;

      const static size_t domainLocalBlockSize = DomainSpaceType::localBlockSize;
      const static size_t rangeLocalBlockSize = RangeSpaceType::localBlockSize;

      typedef Stencil<DomainSpaceType,RangeSpaceType> StencilType;

    private:
      typedef PetscSlaveDofProvider< DomainSpaceType > RowPetscSlaveDofsType;
      typedef PetscSlaveDofProvider< RangeSpaceType  > ColPetscSlaveDofsType;
      enum Status {statAssembled=0,statAdd=1,statInsert=2,statGet=3};

    public:
      typedef typename ColPetscSlaveDofsType :: PetscDofMappingType   ColDofMappingType;
      typedef typename RowPetscSlaveDofsType :: PetscDofMappingType   RowDofMappingType;

      // the local matrix
      class LocalMatrix;

      struct LocalMatrixTraits
      {
        typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
        typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;
        typedef LocalMatrix LocalMatrixType;
        typedef typename RangeSpaceType::RangeFieldType RangeFieldType;

        // copied this typedef from spmatrix.hh
        typedef RangeFieldType LittleBlockType;
      };

      //! type of local matrix object
      typedef LocalMatrix  ObjectType;
      typedef ThisType     LocalMatrixFactoryType;
      typedef ObjectStack< LocalMatrixFactoryType > LocalMatrixStackType;

      //! type of local matrix using stacking mechanism
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;
      typedef ColumnObject< ThisType > LocalColumnObjectType;

      /*
       * ctors, dtor, methods...
       */
      PetscLinearOperator ( const std::string &, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        colSlaveDofs_( rangeSpace_ ),
        rowSlaveDofs_( domainSpace_ ),
        sequence_(-1),
        localMatrixStack_( *this ),
        status_(statAssembled)
      {
      }
      PetscLinearOperator ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        colSlaveDofs_( rangeSpace_ ),
        rowSlaveDofs_( domainSpace_ ),
        sequence_(-1),
        localMatrixStack_( *this ),     
        status_(statAssembled)
      {
      }

      //! destructor deleting PETSc Mat object 
      ~PetscLinearOperator ()
      {
        ::Dune::Petsc::MatDestroy( &petscMatrix_ );
      }

      void communicate ()
      {
        ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FINAL_ASSEMBLY );
        ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FINAL_ASSEMBLY );
         status_ = statAssembled;
      }

      const DomainSpaceType& domainSpace () const { return domainSpace_; }
      const RangeSpaceType& rangeSpace () const { return rangeSpace_; }

      void apply ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        ::Dune::Petsc::MatMult( petscMatrix_, *arg.petscVec() , *dest.petscVec() );
      }

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        apply( arg, dest );
      }

      void reserve (const StencilType &stencil) 
      {
        if(sequence_ != domainSpace().sequence())  
        {
          /*
          * initialize the row and column petsc dof mappings
          */
          const PetscInt localRows =
            rowDofMapping().numOwnedDofBlocks() * domainLocalBlockSize;
          const PetscInt localCols =
            colDofMapping().numOwnedDofBlocks() * rangeLocalBlockSize;

          assert( domainLocalBlockSize == rangeLocalBlockSize );
          // create matrix 
          ::Dune::Petsc::MatCreate( &petscMatrix_ );
        
          if( domainLocalBlockSize > 1 ) 
          {
            ::Dune::Petsc::MatSetType( petscMatrix_, MATBAIJ );
            // set block size 
            PetscInt bs = domainLocalBlockSize ;
            ::Dune::Petsc::MatSetBlockSize( petscMatrix_, bs );
          }
          else 
          {
            ::Dune::Petsc::MatSetType( petscMatrix_, MATAIJ );
          }

          // set sizes of the matrix 
          ::Dune::Petsc::MatSetSizes( petscMatrix_, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE );

          std::vector<int> d_nnz(localRows,0);
          typedef typename StencilType::GlobalStencilType GlobalStencilType;
          typedef typename GlobalStencilType::const_iterator StencilIteratorType;
          const GlobalStencilType &glStencil = stencil.globalStencil();
          StencilIteratorType end = glStencil.end();
          // std::cout << "localRows: " << localRows << std::endl;
          for ( StencilIteratorType it = glStencil.begin(); it != end; ++it)
          {
            int femIndex = it->first;
            int nz = it->second.size();
            int petscIndex = rowDofMapping().localSlaveMapping( femIndex );
              // rowDofMapping().globalMapping( femIndex ) - rowDofMapping().processStartIndex();
            /*
            std::cout << femIndex << " , " 
                      << rowDofMapping().globalMapping( femIndex ) << " , "
                      << petscIndex
                      << " = " << nz
                      << std::endl;
            */
            if ( ! rowDofMapping().isSlave( femIndex ) )
            {
              // std::cout << "inserted..." << std::endl;
              assert( petscIndex >= 0 );
              assert( petscIndex < d_nnz.size() );
              d_nnz[ petscIndex ] = nz;
            }
          }
          ::Dune::Petsc::MatSetUp( petscMatrix_, &d_nnz[0] );
          // ::Dune::Petsc::MatSetUp( petscMatrix_, stencil.maxNZ() );
        } 

        // ::Dune::Petsc::MatSetUp( petscMatrix_ );
        status_ = statAssembled;
      } 

      void clear ()
      {
        ::Dune::Petsc::MatZeroEntries( petscMatrix_ );
      }

      //! interface method from LocalMatrixFactory 
      ObjectType* newObject() const
      {
        return new ObjectType( *this, domainSpace_, rangeSpace_ );
      }
      
      //! return local matrix representation 
      LocalMatrixType localMatrix ( const RowEntityType &rowEntity, const ColumnEntityType &colEntity ) const
      {
        return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
      }
      LocalColumnObjectType localColumn( const ColumnEntityType &colEntity ) const
      {
        return LocalColumnObjectType ( *this, colEntity );
      }

      // just here for debugging
      void view () const 
      {
        ::Dune::Petsc::MatView( petscMatrix_, PETSC_VIEWER_STDOUT_WORLD );
      }
      /* Not tested yet
      void viewMatlab (const char *filename) const
      {
        ::Dune::Petsc::MatViewMatlab( petscMatrix_, filename );
      }
      */

      // return reference to PETSc matrix object 
      Mat& petscMatrix () const { return petscMatrix_; }

      //! return reference to row global mapping 
      const RowDofMappingType& rowDofMapping() const { return rowSlaveDofs_.dofMapping(); }
      //! return reference to column global mapping 
      const ColDofMappingType& colDofMapping() const { return colSlaveDofs_.dofMapping(); }

    private:
      PetscLinearOperator ();

      void setStatus(const Status &newstatus) const
      {
#if 0
        if (status_ != statAssembled && status_ != newstatus)
        {
          ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FLUSH_ASSEMBLY );
          ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FLUSH_ASSEMBLY );
        }
#endif
        status_ = newstatus;
      }

      /*
       * data fields
       */
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;
      ColPetscSlaveDofsType colSlaveDofs_;
      RowPetscSlaveDofsType rowSlaveDofs_;

      int sequence_;
      mutable Mat petscMatrix_;

      mutable LocalMatrixStackType localMatrixStack_;
      mutable Status status_;
    };



    /* ========================================
     * class PetscLinearOperator::LocalMatrix
     */
    template< typename DomainFunction, typename RangeFunction >
    class PetscLinearOperator< DomainFunction, RangeFunction >::LocalMatrix 
    : public LocalMatrixDefault< typename PetscLinearOperator< DomainFunction, RangeFunction >::LocalMatrixTraits >
    {
      typedef LocalMatrix ThisType;
      typedef LocalMatrixDefault< typename PetscLinearOperator< DomainFunction, RangeFunction >::LocalMatrixTraits >  BaseType;

      typedef PetscLinearOperator< DomainFunction, RangeFunction > PetscLinearOperatorType;


    public:
      typedef PetscInt                                            DofIndexType;
      typedef std::vector< DofIndexType >                         IndexVectorType;
      typedef typename DomainFunction::DiscreteFunctionSpaceType  DomainSpaceType;
      typedef typename RangeFunction::DiscreteFunctionSpaceType   RangeSpaceType;
      typedef typename DomainSpaceType::BasisFunctionSetType      DomainBasisFunctionSetType;
      typedef typename RangeSpaceType::BasisFunctionSetType       RangeBasisFunctionSetType;

      enum { littleCols  = RangeSpaceType::localBlockSize };
      enum { littleRows  = DomainSpaceType ::localBlockSize };

    private:

      // needed for .mapEach below
      template< typename PetscMapping >
      struct PetscAssignFunctor
      {
        explicit PetscAssignFunctor ( const PetscMapping &petscMapping, IndexVectorType &indices )
        : petscMapping_( petscMapping ),
          indices_( indices )
        {}

        template< typename T >
        void operator() ( const std::size_t localIndex, T globalIndex ) { indices_[ localIndex ] = petscMapping_.globalMapping( globalIndex ); }

      private:
        const PetscMapping &petscMapping_;
        IndexVectorType &indices_;
      };

    public:

      LocalMatrix ( const PetscLinearOperatorType &petscLinOp, 
                    const DomainSpaceType &domainSpace, 
                    const RangeSpaceType &rangeSpace )
      : BaseType( domainSpace, rangeSpace ),
        petscLinearOperator_( petscLinOp ),
        values_(0)
      {}

      void finalize()
      {
        return; //!
        // USE: MatSetValuesBlocked for bs>1!
        if (status_ == statAdd) // only add is cached at the moment 
        {
          PetscScalar v[rows()][columns()];
          for (int r=0;r<rows();++r)
          {
            for (int c=0;c<columns();++c)
            {
              int idx = c + r*columns();
              v[r][c] = values_[idx];
            }
          }
          ::Dune::Petsc::MatSetValues( petscMatrix(), rows(), &(rowIndices_[0]), columns(), &(colIndices_[0]), 
                                       &v[0][0], ADD_VALUES );
        }
      }

      void init ( const RowEntityType &rowEntity, const ColumnEntityType &colEntity ) 
      {
        // call initialize on base class 
        BaseType :: init( rowEntity, colEntity );

        // setup row indices and also store number of local rows 
        setupIndices( domainSpace().blockMapper(),  petscLinearOperator_.rowDofMapping(), rowEntity, littleRows, rowIndices_ );

        // setup col indices and also store number of local cols 
        setupIndices( rangeSpace().blockMapper(), petscLinearOperator_.colDofMapping(), colEntity, littleCols, colIndices_ );
        
        values_.resize( columns()*rows(), 0. );
        for (int r=0;r<rows();++r)
        {
          for (int c=0;c<columns();++c)
          {
            int idx = c + r*columns();
            values_[idx] = 0;
          }
        }
        status_ = statAssembled;
        petscLinearOperator_.setStatus(status_);
      }

      inline void add ( const int localRow, const int localCol, const RangeFieldType &value )
      {
        assert( status_==statAssembled || status_==statAdd );
        status_ = statAdd;
        petscLinearOperator_.setStatus(status_);
        ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) , value, ADD_VALUES );
        //! int idx = localCol + localRow*columns();
        //! assert( idx < values_.size() );
        //! values_[idx] += value;
      }
      inline void set(const int localRow, const int localCol, const RangeFieldType &value )
      {
        assert( status_==statAssembled || status_==statInsert );
        status_ = statInsert;
        petscLinearOperator_.setStatus(status_);
        ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) , value, INSERT_VALUES );
        // int idx = localCol + localRow*columns();
        // assert( idx < values_.size() );
        // values_[idx] = value;
      }
      //! set matrix row to zero
      void clearRow ( const int localRow )
      {
        assert( status_==statAssembled || status_==statInsert );
        status_ = statInsert;
        petscLinearOperator_.setStatus(status_);
        const int col = this->columns();
        for(int localCol=0; localCol<col; ++localCol) 
          ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) ,
              (localCol==localRow)?1:0., INSERT_VALUES );
        /*
        ::Dune::Petsc::MatAssemblyBegin( petscMatrix(), MAT_FLUSH_ASSEMBLY );
        ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix(), MAT_FLUSH_ASSEMBLY );
        const int r[] = {globalRowIndex( localRow )};
        ::MatZeroRows(petscMatrix(),1,r,0,0,0);
        */
      }
      inline const RangeFieldType get ( const int localRow, const int localCol ) const
      {
        assert( status_==statAssembled || status_==statGet );
        status_ = statGet;
        petscLinearOperator_.setStatus(status_);
        RangeFieldType v[1];
        const int r[] = {globalRowIndex( localRow )};
        const int c[] = {globalColIndex( localCol )};
        ::Dune::Petsc::MatGetValues( petscMatrix(), 1, r, 1, c, v );
        return v[0];
      }

    private:
      LocalMatrix ();

      Mat& petscMatrix () { return petscLinearOperator_.petscMatrix_; }
      const Mat& petscMatrix () const { return petscLinearOperator_.petscMatrix_; }

      // Used to setup row/column indices. DofMapper is the DUNE DoF mapper
      template< typename DofMapper, typename PetscMapping, typename Entity > 
      void setupIndices ( const DofMapper &dofMapper, const PetscMapping &petscMapping, const Entity &entity, 
                          const int blockSize, IndexVectorType &indices )
      {
        const int blockDofs = dofMapper.numDofs( entity ) ;
        const int numDofs   = blockDofs * blockSize ;  
        blockIndices_.resize( blockDofs );
        indices.resize( numDofs );
        // map global dofs (blocked)
        dofMapper.mapEach( entity, PetscAssignFunctor< PetscMapping >( petscMapping, blockIndices_ ) );
        // compute non blocked dofs 
        for( int b=0, dof=0; b<blockDofs; ++ b) 
        {
          int globalDof = blockIndices_[ b ] * blockSize ; 
          for( int d=0; d<blockSize; ++d, ++dof, ++globalDof ) 
          {
            indices[ dof ] = globalDof;
          }
        }
      }

    public:
      const int rows()    const { return rowIndices_.size(); }
      const int columns() const { return colIndices_.size(); }


    private:
      DofIndexType globalRowIndex( const int localRow ) const 
      { 
        assert( localRow < static_cast< int >( rowIndices_.size() ) );
        return rowIndices_[ localRow ]; 
      }

      DofIndexType globalColIndex( const int localCol ) const 
      { 
        assert( localCol < static_cast< int >( colIndices_.size() ) );
        return colIndices_[ localCol ]; 
      }

      /*
       * data fields
       */
      const PetscLinearOperatorType &petscLinearOperator_;
      IndexVectorType blockIndices_;
      IndexVectorType rowIndices_;
      IndexVectorType colIndices_;
      std::vector<RangeFieldType> values_;
      Status status_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #if defined HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
