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

#if defined HAVE_PETSC

#include "petscmat.h"


namespace Dune 
{
  namespace Fem 
  {

    /* ========================================
     * class PetscLinearOperator
     */
    template< typename DomainSpace, typename RangeSpace >
    class PetscLinearOperator 
    {
      typedef PetscLinearOperator< DomainSpace, RangeSpace > ThisType;
    public:

      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;
      typedef typename RangeSpaceType::RangeFieldType RangeFieldType;
      typedef PetscDiscreteFunction< DomainSpaceType > DomainFunctionType;
      typedef PetscDiscreteFunction< RangeSpaceType > RangeFunctionType;
      typedef typename DomainSpaceType::GridPartType::template Codim< 0 >::EntityType ColEntityType;
      typedef typename RangeSpaceType::GridPartType::template Codim< 0 >::EntityType RowEntityType;

      const static size_t domainLocalBlockSize = DomainSpaceType::localBlockSize;
      const static size_t rangeLocalBlockSize = RangeSpaceType::localBlockSize;

    private:
      typedef PetscSlaveDofProvider< DomainSpaceType > ColPetscSlaveDofsType;
      typedef PetscSlaveDofProvider< RangeSpaceType  > RowPetscSlaveDofsType;

    public:
      typedef typename ColPetscSlaveDofsType :: PetscDofMappingType   ColDofMappingType;
      typedef typename RowPetscSlaveDofsType :: PetscDofMappingType   RowDofMappingType;

      // the local matrix
      class LocalMatrix;

      struct LocalMatrixTraits
      {
        typedef DomainSpace DomainSpaceType;
        typedef RangeSpace RangeSpaceType;
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

      /*
       * ctors, dtor, methods...
       */
      PetscLinearOperator ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        colSlaveDofs_( domainSpace_ ),
        rowSlaveDofs_( rangeSpace_ ),
        localMatrixStack_( *this )     
      {
        /*
         * initialize the row and column petsc dof mappings
         */
        const PetscInt localRows = rowDofMapping().numOwnedDofBlocks();
        const PetscInt localCols = colDofMapping().numOwnedDofBlocks();

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
      }

      //! destructor deleting PETSc Mat object 
      ~PetscLinearOperator ()
      {
        ::Dune::Petsc::MatDestroy( &petscMatrix_ );
      }

      // TODO: This is not an interface method. Remove it? But right now I need it.
      void communicate ()
      {
        ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FINAL_ASSEMBLY );
        ::Dune::Petsc::MatAssemblyEnd( petscMatrix_, MAT_FINAL_ASSEMBLY );
      }

      const DomainSpaceType& domainSpace () const { return domainSpace_; }
      const RangeSpaceType& rangeSpace () const { return rangeSpace_; }

      void apply ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        ::Dune::Petsc::MatMult( petscMatrix_, *arg.petscVec() , *dest.petscVec() );
      }

      void reserve () {} 

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
      LocalMatrixType localMatrix ( const RowEntityType &rowEntity, const ColEntityType &colEntity )
      {
        return LocalMatrixType(localMatrixStack_,rowEntity,colEntity);
      }

      // just here for debugging
      #ifndef NDEBUG
      void view () const 
      {
        ::Dune::Petsc::MatView( petscMatrix_, PETSC_VIEWER_STDOUT_WORLD );
      }
      #endif

      // return reference to PETSc matrix object 
      Mat& petscMatrix () const { return petscMatrix_; }

      //! return reference to row global mapping 
      const RowDofMappingType& rowDofMapping() const { return rowSlaveDofs_.dofMapping(); }
      //! return reference to column global mapping 
      const ColDofMappingType& colDofMapping() const { return colSlaveDofs_.dofMapping(); }

    private:
      PetscLinearOperator ();

      /*
       * data fields
       */
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;
      ColPetscSlaveDofsType colSlaveDofs_;
      RowPetscSlaveDofsType rowSlaveDofs_;

      mutable Mat petscMatrix_;

      mutable LocalMatrixStackType localMatrixStack_;
    };



    /* ========================================
     * class PetscLinearOperator::LocalMatrix
     */
    template< typename DomainSpace, typename RangeSpace >
    class PetscLinearOperator< DomainSpace, RangeSpace >::LocalMatrix 
    : public LocalMatrixDefault< typename PetscLinearOperator< DomainSpace, RangeSpace >::LocalMatrixTraits >
    {
      typedef LocalMatrix ThisType;
      typedef LocalMatrixDefault< typename PetscLinearOperator< DomainSpace, RangeSpace >::LocalMatrixTraits >  BaseType;

      typedef PetscLinearOperator< DomainSpace, RangeSpace > PetscLinearOperatorType;


    public:
      typedef PetscInt                                        DofIndexType;
      typedef std::vector< DofIndexType >                     IndexVectorType;
      typedef DomainSpace                                     DomainSpaceType;
      typedef RangeSpace                                      RangeSpaceType;
      typedef typename DomainSpaceType::BaseFunctionSetType   DomainBaseFunctionSetType;
      typedef typename RangeSpaceType::BaseFunctionSetType    RangeBaseFunctionSetType;

      enum { littleCols  = DomainSpaceType::localBlockSize };
      enum { littleRows  = RangeSpaceType ::localBlockSize };

      typedef Dune::FieldMatrix< typename DomainSpaceType :: RangeFieldType, 
                  littleRows, littleCols > LittleBlockType ;

      typedef Dune::DynamicMatrix< LittleBlockType > DynamicMatrixType;

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
        localMatrix_( rangeSpace.blockMapper().maxNumDofs(), 
                      domainSpace.blockMapper().maxNumDofs() )
      {}

      void init ( const RowEntityType &rowEntity, const ColEntityType &colEntity ) 
      {
        // call initialize on base class 
        BaseType :: init( rowEntity, colEntity );

        if( littleRows > 1 )
        {
          // set all values to 0 
          localMatrix_ = 0;
        }

        setupIndices( rangeSpace().blockMapper(), petscLinearOperator_.rowDofMapping(), rowEntity, rowIndices_ );
        setupIndices( domainSpace().blockMapper(), petscLinearOperator_.colDofMapping(), colEntity, colIndices_ );
      }

      inline void add ( const int localRow, const int localCol, const RangeFieldType &value )
      {
        if( littleRows > 1 ) 
        {
          const int row = (int) localRow / littleRows;
          const int col = (int) localCol / littleCols;
          const int lRow = localRow%littleRows;
          const int lCol = localCol%littleCols;

          // add value to local storage 
          localMatrix_[ row ][ col ][ lRow ][ lCol ] += value ;
        }
        else
         ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) , value, ADD_VALUES );
      }

      ~LocalMatrix() 
      {
        if( littleRows > 1 ) 
        {
          const int colSize = columns() * littleCols ;

          PetscScalar v[ colSize ];
          PetscInt c[ colSize ] ;

          const int nrows = rows();
          const int cols = columns();
          for( int row = 0 ; row < nrows; ++ row ) 
          {
            // set row indices 
            PetscInt idxm[ 1 ] = { PetscInt( globalRowIndex( row ) * littleRows ) };
            for(int lr = 0; lr < littleRows; ++lr ) 
            {
              for(int col = 0, colG = 0; col < cols ; ++col ) 
              {
                PetscInt colIdx = globalColIndex( col ) * littleCols ;
                for( int lc = 0; lc < littleCols ; ++ lc, ++ colG, ++colIdx ) 
                {
                  v[ colG ] = localMatrix_[ row ][ col ][ lr ][ lc ];
                  c[ colG ] = colIdx ;
                }
              }

              ::Dune::Petsc::MatSetValues( petscMatrix(), 1, idxm, colSize, c, v, ADD_VALUES );
            }
          }
        }
      }

    private:
      LocalMatrix ();


      Mat& petscMatrix () { return petscLinearOperator_.petscMatrix_; }
      const Mat& petscMatrix () const { return petscLinearOperator_.petscMatrix_; }

      // Used to setup row/column indices. DofMapper is the DUNE DoF mapper
      template< typename DofMapper, typename PetscMapping, typename Entity > 
      void setupIndices ( const DofMapper &dofMapper, const PetscMapping &petscMapping, const Entity &entity, IndexVectorType &indices )
      {
        indices.resize( dofMapper.numDofs( entity ) );
        dofMapper.mapEach( entity, PetscAssignFunctor< PetscMapping >( petscMapping, indices ) );
      }

      const int rows() const { return rowIndices_.size(); }
      const int columns() const { return colIndices_.size(); }

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
      IndexVectorType rowIndices_;
      IndexVectorType colIndices_;

      DynamicMatrixType localMatrix_; 
    };


  } // namespace Fem

} // namespace Dune

#endif // #if defined HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
