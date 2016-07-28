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

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/operator/matrix/functor.hh>

#include <dune/fem/storage/objectstack.hh>

#if HAVE_PETSC

#include "petscmat.h"


namespace Dune
{
  namespace Fem
  {

    struct PetscMatrixParameter
      : public MatrixParameter
    {
      typedef MatrixParameter BaseType;

      PetscMatrixParameter( const std::string keyPrefix = "petscmatrix." )
        : BaseType( keyPrefix )
      {}

    };

    /* ========================================
     * class PetscLinearOperator
     */
    template< typename DomainFunction, typename RangeFunction >
    class PetscLinearOperator
    : public Fem::Operator< DomainFunction, RangeFunction >
    {
      typedef PetscLinearOperator< DomainFunction, RangeFunction > ThisType;
    public:
      typedef Mat MatrixType;
      typedef DomainFunction DomainFunctionType;
      typedef RangeFunction RangeFunctionType;
      typedef typename DomainFunctionType::RangeFieldType DomainFieldType;
      typedef typename RangeFunctionType::RangeFieldType RangeFieldType;
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

      typedef PetscDiscreteFunction< DomainSpaceType > PetscDomainFunctionType;
      typedef PetscDiscreteFunction< RangeSpaceType  > PetscRangeFunctionType;

      typedef typename DomainSpaceType::GridPartType::template Codim< 0 >::EntityType RowEntityType;
      typedef typename RangeSpaceType::GridPartType::template Codim< 0 >::EntityType ColumnEntityType;

      const static size_t domainLocalBlockSize = DomainSpaceType::localBlockSize;
      const static size_t rangeLocalBlockSize = RangeSpaceType::localBlockSize;

      static_assert( domainLocalBlockSize == rangeLocalBlockSize, "PetscLinearOperator only works for domainLocalBlockSize == rangeLocalBlockSize. " );

    private:
      typedef PetscSlaveDofProvider< DomainSpaceType > RowPetscSlaveDofsType;
      typedef PetscSlaveDofProvider< RangeSpaceType  > ColPetscSlaveDofsType;
      enum Status {statAssembled=0,statAdd=1,statInsert=2,statGet=3,statNothing=4};

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
      PetscLinearOperator ( const std::string &, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace,
                            const MatrixParameter& param = PetscMatrixParameter() )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        colSlaveDofs_( rangeSpace_ ),
        rowSlaveDofs_( domainSpace_ ),
        sequence_(-1),
        localMatrixStack_( *this ),
        status_(statNothing)
      {
      }
      PetscLinearOperator ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace,
                            const MatrixParameter& param = PetscMatrixParameter() )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        colSlaveDofs_( rangeSpace_ ),
        rowSlaveDofs_( domainSpace_ ),
        sequence_(-1),
        localMatrixStack_( *this ),
        status_(statNothing)
      {
      }

      //! destructor deleting PETSc Mat object
      ~PetscLinearOperator ()
      {
        removeObj();
      }

      void communicate ()
      {
        ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FINAL_ASSEMBLY );
        ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FINAL_ASSEMBLY );
         status_ = statAssembled;
      }

      const DomainSpaceType& domainSpace () const { return domainSpace_; }
      const RangeSpaceType& rangeSpace () const { return rangeSpace_; }

      /** \brief application operator for arbitrary DiscreteFunction
       *  \note This functions needs to make copies of the dof vectors into
       *  PetscDiscreteFunction */
      template <class DF, class RF>
      void apply ( const DF &arg, RF &dest ) const
      {
        if( ! petscArg_ )
          petscArg_.reset( new PetscDomainFunctionType( "PetscOp-arg", domainSpace_ ) );
        if( ! petscDest_ )
          petscDest_.reset( new PetscRangeFunctionType( "PetscOp-arg", rangeSpace_ ) );

        petscArg_->assign( arg );
        ::Dune::Petsc::MatMult( petscMatrix_, *(petscArg_->petscVec()) , *(petscDest_->petscVec()) );
        dest.assign( *petscDest_ );
      }

      /** \brief application operator for PetscDiscreteFunction */
      void apply ( const PetscDomainFunctionType &arg, PetscRangeFunctionType &dest ) const
      {
        ::Dune::Petsc::MatMult( petscMatrix_, *arg.petscVec() , *dest.petscVec() );
      }

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        apply( arg, dest );
      }

      void reserve ()
      {
        reserve( SimpleStencil<DomainSpaceType,RangeSpaceType>(0) );
      }

      //! reserve memory for assemble based on the provided stencil
      template <class StencilType>
      void reserve (const StencilType &stencil)
      {
        if(sequence_ != domainSpace().sequence())
        {
          // clear Petsc Mat
          removeObj();

          // reset temporary Petsc discrete functions
          petscArg_.reset();
          petscDest_.reset();

          // update dof mappings
          rowSlaveDofs_.update();
          colSlaveDofs_.update();

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

          PetscInt bs = 1;
          if( domainLocalBlockSize > 1 )
          {
            bs = domainLocalBlockSize ;
            ::Dune::Petsc::MatSetType( petscMatrix_, MATBAIJ );
            // set block size
            ::Dune::Petsc::MatSetBlockSize( petscMatrix_, bs );
          }
          else
          {
            ::Dune::Petsc::MatSetType( petscMatrix_, MATAIJ );
          }

          // set sizes of the matrix
          ::Dune::Petsc::MatSetSizes( petscMatrix_, localRows, localCols, PETSC_DETERMINE, PETSC_DETERMINE );

          if (std::is_same< StencilType,SimpleStencil<DomainSpaceType,RangeSpaceType> >::value)
            ::Dune::Petsc::MatSetUp( petscMatrix_, bs, stencil.maxNonZerosEstimate() );
          else
          {
            std::vector<int> d_nnz(localRows/bs,0);
            std::vector<int> o_nnz(localRows/bs,0);
            typedef typename StencilType::GlobalStencilType GlobalStencilType;
            typedef typename GlobalStencilType::const_iterator StencilIteratorType;
            const GlobalStencilType &glStencil = stencil.globalStencil();
            StencilIteratorType end = glStencil.end();
            for ( StencilIteratorType it = glStencil.begin(); it != end; ++it)
            {
              int femIndex = it->first;
              if ( rowDofMapping().isSlave( femIndex ) ) continue;
              // Remark: ghost entities should not be inserted into the stencil for dg to
              // get optimal results but they are needed for istl....
              int nzDiag = 0;
              int nzOff  = 0;

              typedef typename StencilType::LocalStencilType LocalStencilType;
              typedef typename LocalStencilType::const_iterator LocalStencilIteratorType;
              LocalStencilIteratorType endLocal = it->second.end();
              for ( LocalStencilIteratorType itLocal = it->second.begin(); itLocal != endLocal; ++itLocal)
              {
                if (!rowDofMapping().isSlave( *itLocal ))
                {
                  ++nzDiag;
                  // std::cout << "diag: (" << rowDofMapping().localSlaveMapping( femIndex )
                  //   << "," << rowDofMapping().localSlaveMapping( *itLocal )
                  //   << ") ";
                }
                else
                  ++nzOff;
              }
              // std::cout << std::endl;
              // std::cout << "nz: " << rowDofMapping().localSlaveMapping( femIndex )
              //           << " " << nzDiag << " " << nzOff << std::endl;

              int petscIndex = rowDofMapping().localSlaveMapping( femIndex );
              assert( petscIndex >= 0 );
              assert( petscIndex < ( int ) d_nnz.size() );
              d_nnz[petscIndex] = nzDiag;
              o_nnz[petscIndex] = nzOff;
            }
            ::Dune::Petsc::MatSetUp( petscMatrix_, bs, &d_nnz[0], &o_nnz[0] );
          }
          sequence_ = domainSpace().sequence();
        }

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

      template< class LocalMatrix >
      void addLocalMatrix ( const RowEntityType &domainEntity, const ColumnEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        assert( status_==statAssembled || status_==statAdd );
        setStatus( statAdd );

        std::vector< PetscInt > r, c;
        setupIndices( rangeSpace_.blockMapper(), rowDofMapping(), rangeEntity, rangeLocalBlockSize, r);
        setupIndices( domainSpace_.blockMapper(), colDofMapping(), domainEntity, domainLocalBlockSize, c);

        std::vector< PetscScalar > v( r.size() * c.size() );
        for( std::size_t i =0 ; i< r.size(); ++i )
          for( std::size_t j =0; j< c.size(); ++j )
            v[ i * c.size() +j ] = localMat.get( i, j );

        ::Dune::Petsc::MatSetValues( petscMatrix_, r.size(), r.data(), c.size(), c.data(), v.data(), ADD_VALUES );
        setStatus( statAssembled );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const RowEntityType &domainEntity, const ColumnEntityType &rangeEntity, const LocalMatrix &localMat, const Scalar &s )
      {
        assert( status_==statAssembled || status_==statAdd );
        setStatus( statAdd );

        std::vector< PetscInt > r, c;
        setupIndices( rangeSpace_.blockMapper(), rowDofMapping(), rangeEntity, rangeLocalBlockSize, r);
        setupIndices( domainSpace_.blockMapper(), colDofMapping(), domainEntity, domainLocalBlockSize, c);

        std::vector< PetscScalar > v( r.size() * c.size() );
        for( std::size_t i =0 ; i< r.size(); ++i )
          for( std::size_t j =0; j< c.size(); ++j )
            v[ i * c.size() +j ] = s * localMat.get( i, j );

        ::Dune::Petsc::MatSetValues( petscMatrix_, r.size(), r.data(), c.size(), c.data(), v.data(), ADD_VALUES );
        setStatus( statAssembled );
      }

      template< class LocalMatrix >
      void setLocalMatrix ( const RowEntityType &domainEntity, const ColumnEntityType &rangeEntity, const LocalMatrix &localMat )
      {
        assert( status_==statAssembled || status_==statInsert );
        setStatus( statInsert );

        std::vector< PetscInt > r, c;
        setupIndices( rangeSpace_.blockMapper(), rowDofMapping(), rangeEntity, rangeLocalBlockSize, r);
        setupIndices( domainSpace_.blockMapper(), colDofMapping(), domainEntity, domainLocalBlockSize, c);

        std::vector< PetscScalar > v( r.size() * c.size() );
        for( std::size_t i =0 ; i< r.size(); ++i )
          for( std::size_t j =0; j< c.size(); ++j )
            v[ i * c.size() +j ] = localMat.get( i, j );

        ::Dune::Petsc::MatSetValues( petscMatrix_, r.size(), r.data(), c.size(), c.data(), v.data(), INSERT_VALUES );

        setStatus( statAssembled );
      }


      template< class LocalMatrix >
      void getLocalMatrix ( const RowEntityType &domainEntity, const ColumnEntityType &rangeEntity, LocalMatrix &localMat ) const
      {
        assert( status_==statAssembled || status_==statGet );
        setStatus( statGet );

        std::vector< PetscInt > r, c;
        setupIndices( rangeSpace_.blockMapper(), rowDofMapping(), rangeEntity, rangeLocalBlockSize, r);
        setupIndices( domainSpace_.blockMapper(), colDofMapping(), domainEntity, domainLocalBlockSize, c);

        std::vector< PetscScalar > v( r.size() * c.size() );
        ::Dune::Petsc::MatGetValues( petscMatrix_, r.size(), r.data(), c.size(), c.data(), v.data() );

        for( std::size_t i =0 ; i< r.size(); ++i )
          for( std::size_t j =0; j< c.size(); ++j )
            localMat.set( i, j, v[ i * c.size() +j ] );

        setStatus( statAssembled );
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

      //! destructor deleting PETSc Mat object
      void removeObj ()
      {
        if( status_ != statNothing )
          ::Dune::Petsc::MatDestroy( &petscMatrix_ );
      }

      void setStatus(const Status &newstatus) const
      {
        status_ = newstatus;
      }

      // Used to setup row/column indices. DofMapper is the DUNE DoF mapper
      template< typename DofMapper, typename PetscMapping, typename Entity >
      static void setupIndices ( const DofMapper &dofMapper, const PetscMapping &petscMapping, const Entity &entity,
          PetscInt blockSize, std::vector< PetscInt > &indices )
      {
        const int blockDofs = dofMapper.numDofs( entity );
        const int numDofs   = blockDofs * blockSize;

        indices.resize( numDofs );

        auto functor = [ &petscMapping, blockSize, &indices  ] ( PetscInt local, PetscInt global )
        {
          const PetscInt block = petscMapping.globalMapping( global );
          for( PetscInt b = 0; b < blockSize; ++b )
            indices[ local *blockSize + b ] = block * blockSize + b;

        };

        // map global dofs (blocked)
        dofMapper.mapEach( entity, functor );
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

      mutable std::unique_ptr< PetscDomainFunctionType > petscArg_;
      mutable std::unique_ptr< PetscRangeFunctionType  > petscDest_;
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

      enum { rangeBlockSize = RangeSpaceType::localBlockSize };
      enum { domainBlockSize = DomainSpaceType ::localBlockSize };

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

      void init ( const RowEntityType &domainEntity, const ColumnEntityType &rangeEntity )
      {
        // call initialize on base class
        BaseType :: init( domainEntity, rangeEntity );

        //*************************************************
        //  The rows belong to the domain space
        //  it's indices are determained by the rangeSpace
        //
        //  The columns belong to the range space
        //  it's indices are determained by the domainSpace
        //*************************************************

        // setup row indices and also store number of local rows
        setupIndices( rangeSpace().blockMapper(), petscLinearOperator_.colDofMapping(), rangeEntity, rangeBlockSize, rowIndices_ );

        // setup col indices and also store number of local cols
        setupIndices( domainSpace().blockMapper(), petscLinearOperator_.rowDofMapping(), domainEntity, domainBlockSize, colIndices_ );


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
      IndexVectorType rowIndices_;
      IndexVectorType colIndices_;
      std::vector<RangeFieldType> values_;
      mutable Status status_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
