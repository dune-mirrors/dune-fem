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

      typedef typename DomainSpaceType::GridPartType::template Codim< 0 >::EntityType ColumnEntityType;
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
        colSlaveDofs_( domainSpace_ ),
        rowSlaveDofs_( rangeSpace_ ),
        sequence_(-1),
        localMatrixStack_( *this )     
      {
      }
      PetscLinearOperator ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        colSlaveDofs_( domainSpace_ ),
        rowSlaveDofs_( rangeSpace_ ),
        sequence_(-1),
        localMatrixStack_( *this )     
      {
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
        ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FINAL_ASSEMBLY );
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

      void reserve () 
      {
        if(sequence_ != domainSpace().sequence())  
        {
          /*
          * initialize the row and column petsc dof mappings
          */
          const PetscInt localRows =
            rowDofMapping().numOwnedDofBlocks() * rangeLocalBlockSize;
          const PetscInt localCols =
            colDofMapping().numOwnedDofBlocks() * domainLocalBlockSize;

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

        ::Dune::Petsc::MatSetUp( petscMatrix_ );
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

      int sequence_;
      mutable Mat petscMatrix_;

      mutable LocalMatrixStackType localMatrixStack_;
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

      enum { littleCols  = DomainSpaceType::localBlockSize };
      enum { littleRows  = RangeSpaceType ::localBlockSize };

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
        petscLinearOperator_( petscLinOp )
      {}

      void init ( const RowEntityType &rowEntity, const ColumnEntityType &colEntity ) 
      {
        // call initialize on base class 
        BaseType :: init( rowEntity, colEntity );

        // setup row indices and also store number of local rows 
        setupIndices( rangeSpace().blockMapper(),  petscLinearOperator_.rowDofMapping(), rowEntity, littleRows, rowIndices_ );

        // setup col indices and also store number of local cols 
        setupIndices( domainSpace().blockMapper(), petscLinearOperator_.colDofMapping(), colEntity, littleCols, colIndices_ );
      }

      inline void add ( const int localRow, const int localCol, const RangeFieldType &value )
      {
        ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) , value, ADD_VALUES );
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
    };


  } // namespace Fem

} // namespace Dune

#endif // #if defined HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
