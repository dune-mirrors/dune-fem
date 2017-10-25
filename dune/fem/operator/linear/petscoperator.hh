#ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
#define DUNE_FEM_PETSCLINEAROPERATOR_HH

#include <cassert>
#include <cstddef>

#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/fem/function/petscdiscretefunction/petscdiscretefunction.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/operator/matrix/functor.hh>
#include <dune/fem/operator/matrix/spmatrix.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/mapper/petsc.hh>
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
    : public Fem::AssembledOperator< DomainFunction, RangeFunction >
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
      enum Status {statAssembled=0,statAdd=1,statInsert=2,statGet=3,statNothing=4};

      typedef PetscMappers< DomainSpaceType > DomainMappersType;
      typedef PetscMappers< RangeSpaceType > RangeMappersType;

      typedef SlaveDofs< typename DomainSpaceType::GridPartType, typename DomainMappersType::GhostMapperType > DomainSlaveDofsType;
      typedef SlaveDofs< typename RangeSpaceType::GridPartType, typename RangeMappersType::GhostMapperType > RangeSlaveDofsType;

    public:
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

      PetscLinearOperator ( const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace,
                            const MatrixParameter& param = PetscMatrixParameter() )
        : domainMappers_( domainSpace ),
          rangeMappers_( rangeSpace ),
          localMatrixStack_( *this ),
          status_(statNothing)
      {}

      PetscLinearOperator ( const std::string &, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace,
                            const MatrixParameter& param = PetscMatrixParameter() )
        : PetscLinearOperator( domainSpace, rangeSpace, param )
      {}

      //! destructor deleting PETSc Mat object
      ~PetscLinearOperator ()
      {
        destroy();
      }

      void communicate ()
      {
        ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FINAL_ASSEMBLY );
        ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FINAL_ASSEMBLY );
         status_ = statAssembled;
      }

      const DomainSpaceType &domainSpace () const { return domainMappers_.space(); }
      const RangeSpaceType &rangeSpace () const { return rangeMappers_.space(); }

      /** \brief application operator for arbitrary DiscreteFunction
       *  \note This functions needs to make copies of the dof vectors into
       *  PetscDiscreteFunction */
      template <class DF, class RF>
      void apply ( const DF &arg, RF &dest ) const
      {
        if( ! petscArg_ )
          petscArg_.reset( new PetscDomainFunctionType( "PetscOp-arg", domainSpace() ) );
        if( ! petscDest_ )
          petscDest_.reset( new PetscRangeFunctionType( "PetscOp-arg", rangeSpace() ) );

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
        domainMappers_.update();
        rangeMappers_.update();

        if(sequence_ != domainSpace().sequence())
        {
          // clear Petsc Mat
          destroy();

          // reset temporary Petsc discrete functions
          petscArg_.reset();
          petscDest_.reset();

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
          const PetscInt localRows = domainMappers_.ghostMapper().interiorSize() * domainLocalBlockSize;
          const PetscInt localCols = rangeMappers_.ghostMapper().interiorSize() * rangeLocalBlockSize;
          const PetscInt globalRows = domainMappers_.parallelMapper().size() * domainLocalBlockSize;
          const PetscInt globalCols = rangeMappers_.parallelMapper().size() * rangeLocalBlockSize;
          ::Dune::Petsc::MatSetSizes( petscMatrix_, localRows, localCols, globalRows, globalCols );

          DomainSlaveDofsType domainSlaveDofs( domainMappers_.ghostMapper() );
          RangeSlaveDofsType rangeSlaveDofs( rangeMappers_.ghostMapper() );

          if (std::is_same< StencilType,SimpleStencil<DomainSpaceType,RangeSpaceType> >::value)
            ::Dune::Petsc::MatSetUp( petscMatrix_, bs, stencil.maxNonZerosEstimate() );
          else
          {
            std::vector< PetscInt > d_nnz( localRows / bs, 0 );
            std::vector< PetscInt > o_nnz( localRows / bs, 0 );
            for( const auto entry : stencil.globalStencil() )
            {
              const int petscIndex = rangeMappers_.ghostIndex( entry.first );
              if( rangeSlaveDofs.isSlave( petscIndex ) )
                continue;

              // Remark: ghost entities should not be inserted into the stencil for dg to
              // get optimal results but they are needed for istl....
              d_nnz[ petscIndex ] = o_nnz[ petscIndex ] = 0;
              for( const auto local : entry.second )
              {
                if( !domainSlaveDofs.isSlave( domainMappers_.ghostIndex( local ) ) )
                  ++d_nnz[ petscIndex ];
                else
                  ++o_nnz[ petscIndex ];
              }
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

      template <class Vector>
      void setUnitRows( const Vector &rows )
      {
        std::vector< PetscInt > r( rows.size() );
        for( std::size_t i =0 ; i< rows.size(); ++i )
        {
          const PetscInt block = rangeMappers_.parallelIndex( rows[ i ] / rangeLocalBlockSize );
          r[ i ] = block * rangeLocalBlockSize + (rows[ i ] % rangeLocalBlockSize);
        }
        ::Dune::Petsc::MatZeroRows( petscMatrix_, r.size(), r.data(), 1. );
      }

      void maskRows ( const RangeFunctionType &maskFunction, PetscScalar diagonal = PetscScalar( 0 ) )
      {
        PetscScalar *maskArray;
        VecGetArray( *maskFunction.getGhostedVector(), &maskArray );

        std::vector< PetscInt > r;
        const PetscInt size = rangeMappers_.parallelMapper().size() * rangeLocalBlockSize;
        for( PetscInt i = 0; i != size; ++i )
        {
          if( maskArray[ i ] != 0 )
            r.push_back( i );
        }
        ::Dune::Petsc::MatZeroRows( petscMatrix_, r.size(), r.data(), diagonal );
      }

      //! interface method from LocalMatrixFactory
      ObjectType* newObject() const
      {
        return new ObjectType( *this, domainSpace(), rangeSpace() );
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
        setupIndices( rangeMappers_, rangeEntity, r );
        setupIndices( domainMappers_, domainEntity, c );

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
        setupIndices( rangeMappers_, rangeEntity, r );
        setupIndices( domainMappers_, domainEntity, c );

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
        setupIndices( rangeMappers_, rangeEntity, r );
        setupIndices( domainMappers_, domainEntity, c );

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
        setupIndices( rangeMappers_, rangeEntity, r );
        setupIndices( domainMappers_, domainEntity, c );

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

      // print matrix just here for debugging
      void print( std::ostream& s ) const
      {
        if( &s == &std::cout || &s == &std::cerr )
        {
          view();
        }
      }

      /* Not tested yet
      void viewMatlab (const char *filename) const
      {
        ::Dune::Petsc::MatViewMatlab( petscMatrix_, filename );
      }
      */

      // return reference to PETSc matrix object
      Mat& petscMatrix () const { return petscMatrix_; }

    private:
      PetscLinearOperator ();

      //! destructor deleting PETSc Mat object
      void destroy ()
      {
        if( status_ != statNothing )
        {
          ::Dune::Petsc::MatDestroy( &petscMatrix_ );
          setStatus( statNothing );
        }
        sequence_ = -1;
      }

      void setStatus(const Status &newstatus) const
      {
        status_ = newstatus;
      }

      template< class DFS, class Entity >
      static void setupIndices ( const PetscMappers< DFS > &mappers, const Entity &entity, std::vector< PetscInt > &indices )
      {
        NonBlockMapper< const typename PetscMappers< DFS >::ParallelMapperType, DFS::localBlockSize > nonBlockMapper( mappers.parallelMapper() );
        nonBlockMapper.map( entity, indices );
      }

      /*
       * data fields
       */
      DomainMappersType domainMappers_;
      RangeMappersType rangeMappers_;

      int sequence_ = -1;
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

        setupIndices( petscLinearOperator_.rangeMappers_, rangeEntity, rowIndices_ );
        setupIndices( petscLinearOperator_.domainMappers_, domainEntity, colIndices_ );

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
        const int globalRowIdx = globalRowIndex( localRow );
        for(int localCol=0; localCol<col; ++localCol)
        {
          ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIdx, globalColIndex( localCol ), 0.0, INSERT_VALUES );
        }

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
