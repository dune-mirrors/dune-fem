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
#include <dune/fem/misc/fmatrixconverter.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>
#include <dune/fem/operator/common/localcontribution.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/operator/matrix/functor.hh>
#include <dune/fem/space/mapper/nonblockmapper.hh>
#include <dune/fem/space/mapper/petsc.hh>
#include <dune/fem/storage/objectstack.hh>

#include <dune/fem/solver/parameter.hh>

#if HAVE_PETSC

#include "petscmat.h"


namespace Dune
{
  namespace Fem
  {

    struct PetscSolverParameter : public LocalParameter< SolverParameter, PetscSolverParameter >
    {
      typedef LocalParameter< SolverParameter, PetscSolverParameter >  BaseType;

    public:
      using BaseType :: parameter ;
      using BaseType :: keyPrefix ;

      PetscSolverParameter( const ParameterReader &parameter = Parameter::container() )
        : BaseType( parameter )
      {}

      PetscSolverParameter( const SolverParameter& sp )
        : PetscSolverParameter( sp.keyPrefix(), sp.parameter() )
      {}

      PetscSolverParameter( const std::string &keyPrefix, const ParameterReader &parameter = Parameter::container() )
        : BaseType( keyPrefix, parameter )
      {}

      bool isPetscSolverParameter() const { return true; }

      static const int boomeramg = 0;
      static const int parasails = 1;
      static const int pilut     = 2;

      int hypreMethod() const
      {
        const std::string hyprePCNames[] = { "boomer-amg", "parasails", "pilu-t" };
        int hypreType = 0;
        if (parameter().exists("petsc.preconditioning.method"))
        {
          hypreType = parameter().getEnum( "petsc.preconditioning.hypre.method", hyprePCNames, 0 );
          std::cout << "WARNING: using deprecated parameter 'petsc.preconditioning.hypre.method' use "
              << keyPrefix() << "preconditioning.hypre.method instead\n";
        }
        else
          hypreType = parameter().getEnum( keyPrefix()+"hypre.method", hyprePCNames, 0 );
        return hypreType;
      }

      int superluMethod() const
      {
        const std::string factorizationNames[] = { "petsc", "superlu", "mumps" };
        int factorType = 0;
        if (parameter().exists("petsc.preconditioning.lu.method"))
        {
          factorType = parameter().getEnum( "petsc.preconditioning.lu.method", factorizationNames, 0 );
          std::cout << "WARNING: using deprecated parameter 'petsc.preconditioning.lu.method' use "
              << keyPrefix() << "preconditioning.lu.method instead\n";
        }
        else
          factorType = parameter().getEnum( keyPrefix()+"preconditioning.lu.method", factorizationNames, 0 );
        return factorType;
      }

      bool viennaCL () const {
        return parameter().getValue< bool > ( keyPrefix() + "petsc.viennacl", false );
      }
      bool blockedMode () const {
        return parameter().getValue< bool > ( keyPrefix() + "petsc.blockedmode", true );
      }
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

      typedef typename DomainSpaceType::GridPartType::template Codim< 0 >::EntityType DomainEntityType;
      typedef typename RangeSpaceType::GridPartType::template Codim< 0 >::EntityType  RangeEntityType;

      static const unsigned int domainLocalBlockSize = DomainSpaceType::localBlockSize;
      static const unsigned int rangeLocalBlockSize  = RangeSpaceType::localBlockSize;

      static constexpr bool squareBlocked = domainLocalBlockSize == rangeLocalBlockSize ;
      static constexpr bool blockedMatrix = domainLocalBlockSize > 1 && squareBlocked;

      // is this right? It should be rangeBS x domainBS - the system is
      // Ad=r so domain gives columns and r gives range
      typedef FlatFieldMatrix< RangeFieldType, domainLocalBlockSize, rangeLocalBlockSize > MatrixBlockType;
      typedef MatrixBlockType  block_type;

    private:
      enum Status {statAssembled=0,statAdd=1,statInsert=2,statGet=3,statNothing=4};

      typedef PetscMappers< DomainSpaceType > DomainMappersType;
      typedef PetscMappers< RangeSpaceType > RangeMappersType;

      typedef AuxiliaryDofs< typename DomainSpaceType::GridPartType, typename DomainMappersType::GhostMapperType > DomainAuxiliaryDofsType;
      typedef AuxiliaryDofs< typename RangeSpaceType::GridPartType, typename RangeMappersType::GhostMapperType > RangeAuxiliaryDofsType;

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
                            const PetscSolverParameter& param = PetscSolverParameter() )
        : domainMappers_( domainSpace ),
          rangeMappers_( rangeSpace ),
          localMatrixStack_( *this ),
          status_(statNothing),
          viennaCL_( param.viennaCL() ),
          blockedMode_( blockedMatrix && (!viennaCL_) && param.blockedMode() )
      {}

      PetscLinearOperator ( const std::string &, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace,
                            const PetscSolverParameter& param = PetscSolverParameter() )
        : PetscLinearOperator( domainSpace, rangeSpace, param )
      {}

      //! destructor deleting PETSc Mat object
      ~PetscLinearOperator ()
      {
        destroy();
      }

      void flushAssembly()
      {
        ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FLUSH_ASSEMBLY );
        ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FLUSH_ASSEMBLY );
        // set variable directly since setStatus might be disabled
        status_ = statAssembled;
      }

      void finalize ()
      {
        if( ! finalizedAlready() )
        {
          ::Dune::Petsc::MatAssemblyBegin( petscMatrix_, MAT_FINAL_ASSEMBLY );
          ::Dune::Petsc::MatAssemblyEnd  ( petscMatrix_, MAT_FINAL_ASSEMBLY );

          setUnitRowsImpl( unitRows_, auxRows_ );
          // remove all cached unit rows
          unitRows_.clear();
          auxRows_.clear();
        }
      }

    protected:
      bool finalizedAlready() const
      {
        PetscBool assembled = PETSC_FALSE ;
        ::Dune::Petsc::MatAssembled( petscMatrix_, &assembled );
        return assembled == PETSC_TRUE;
      }

      void finalizeAssembly () const
      {
        const_cast< ThisType& > (*this).finalize();
      }

    public:
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
        apply( *petscArg_, *petscDest_ );
        dest.assign( *petscDest_ );
      }

      /** \brief application operator for PetscDiscreteFunction */
      void apply ( const PetscDomainFunctionType &arg, PetscRangeFunctionType &dest ) const
      {
        // make sure matrix is in correct state
        finalizeAssembly();
        ::Dune::Petsc::MatMult( petscMatrix_, *arg.petscVec() , *dest.petscVec() );
      }

      void operator() ( const DomainFunctionType &arg, RangeFunctionType &dest ) const
      {
        apply( arg, dest );
      }

      void reserve ()
      {
        reserve( SimpleStencil<DomainSpaceType,RangeSpaceType>(0), true );
      }

      template <class Set>
      void reserve (const std::vector< Set >& sparsityPattern )
      {
        reserve( StencilWrapper< DomainSpaceType,RangeSpaceType, Set >( sparsityPattern ), true );
      }

      //! reserve memory for assemble based on the provided stencil
      template <class StencilType>
      void reserve (const StencilType &stencil, const bool isSimpleStencil = false )
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

          // create matrix
          ::Dune::Petsc::MatCreate( domainSpace().gridPart().comm(), &petscMatrix_ );

          // set sizes of the matrix (columns == domain and rows == range)
          const PetscInt localCols = domainMappers_.ghostMapper().interiorSize() * domainLocalBlockSize;
          const PetscInt localRows = rangeMappers_.ghostMapper().interiorSize() * rangeLocalBlockSize;

          const PetscInt globalCols = domainMappers_.parallelMapper().size() * domainLocalBlockSize;
          const PetscInt globalRows = rangeMappers_.parallelMapper().size() * rangeLocalBlockSize;

          ::Dune::Petsc::MatSetSizes( petscMatrix_, localRows, localCols, globalRows, globalCols );

          PetscInt bs = 1;
          if( viennaCL_ )
          {
            ::Dune::Petsc::MatSetType( petscMatrix_, MATAIJVIENNACL );
          }
          else if( blockedMode_ )
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
          /*
          std::cout << "Matrix dimension with bs=" << bs << "   "
            << localRows << "x" << localCols << "   "
            << globalRows << "x" <<  globalCols << "   "
            << rangeLocalBlockSize/bs << "x" << domainLocalBlockSize/bs << "    "
            << std::endl;
          */

          // there is an issue with the stencil and non zero pattern in the
          // case of domainSpace != rangeSpace. In this case additional
          // mallocs are required during matrix assembly which heavily
          // impacts the preformance. As a quick fix we use a global value
          // for the number of non zeros per row. That leads to a
          // significant increase in memory usage but not to much
          // computational overhead in assembly. The issue with the stencil
          // is a bug and needs fixing....
          if( isSimpleStencil || std::is_same< StencilType,SimpleStencil<DomainSpaceType,RangeSpaceType> >::value )
          {
            ::Dune::Petsc::MatSetUp( petscMatrix_, bs, domainLocalBlockSize * stencil.maxNonZerosEstimate() );
          }
          else
          {
            DomainAuxiliaryDofsType domainAuxiliaryDofs( domainMappers_.ghostMapper() );
            RangeAuxiliaryDofsType rangeAuxiliaryDofs( rangeMappers_.ghostMapper() );

            std::vector< PetscInt > d_nnz( localRows / bs, 0 );
            std::vector< PetscInt > o_nnz( localRows / bs, 0 );
            for( const auto& entry : stencil.globalStencil() )
            {
              const int petscIndex = rangeMappers_.ghostIndex( entry.first );
              if( rangeAuxiliaryDofs.contains( petscIndex ) )
                continue;

              for (unsigned int rb = 0; rb<rangeLocalBlockSize/bs; ++rb)
              {
                const size_t row = petscIndex*rangeLocalBlockSize/bs + rb;
                // Remark: ghost entities should not be inserted into the stencil for dg to
                // get optimal results but they are needed for istl....
                assert( row < size_t(localRows/bs) );
                d_nnz[ row ] = o_nnz[ row ] = 0;
                for( const auto local : entry.second )
                {
                  if( !domainAuxiliaryDofs.contains( domainMappers_.ghostIndex( local ) ) )
                    d_nnz[ row ] += domainLocalBlockSize/bs;
                  else
                    o_nnz[ row ] += domainLocalBlockSize/bs;
                }
              }
            }
            ::Dune::Petsc::MatSetUp( petscMatrix_, bs, &d_nnz[0], &o_nnz[0] );
          }
          sequence_ = domainSpace().sequence();
        }

        flushAssembly();
        status_ = statAssembled ;
      }

      void clear ()
      {
        flushAssembly();
        ::Dune::Petsc::MatZeroEntries( petscMatrix_ );
        flushAssembly();

        unitRows_.clear();
        auxRows_.clear();
      }

      void setUnitRowsImpl( const std::vector< PetscInt >& r,
                            const std::vector< PetscInt >& a )
      {
        ::Dune::Petsc::MatZeroRows( petscMatrix_, r.size(), r.data(), 1.0 );
        ::Dune::Petsc::MatZeroRows( petscMatrix_, a.size(), a.data(), 0.0 );
      }

      template <class Container> // could bet std::set or std::vector
      void setUnitRows( const Container& unitRows, const Container& auxRows )
      {
        std::vector< PetscInt > r, a;
        r.reserve( unitRows.size() );
        a.reserve( auxRows.size() );

        auto setupRows = [this] (const Container& rows, std::vector< PetscInt >& r )
        {
          for( const auto& row : rows )
          {
            const PetscInt block = this->rangeMappers_.parallelIndex( row / rangeLocalBlockSize );
            r.push_back( block * rangeLocalBlockSize + (row % rangeLocalBlockSize) );
          }
        };

        setupRows( unitRows, r );
        setupRows( auxRows,  a );

        if( finalizedAlready() )
        {
          setUnitRowsImpl( r, a );
        }
        else
        {
          unitRows_.reserve( unitRows_.size() + r.size() );
          for( const auto& row : r )
            unitRows_.push_back( row );

          auxRows_.reserve( auxRows_.size() + a.size() );
          for( const auto& row : a )
            auxRows_.push_back( row );
        }
      }

      //! interface method from LocalMatrixFactory
      ObjectType* newObject() const
      {
        return new ObjectType( *this, domainSpace(), rangeSpace() );
      }

      /** \deprecated Use TemporaryLocalMatrix in combination with
        *             {add,set,get}LocalMatrix on matrix object
        *  return local matrix object
        */
      [[deprecated("Use TemporaryLocalMatrix,LocalContribution and {get,add,set}LocalMatrix")]]
      LocalMatrixType localMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity ) const
      {
        return LocalMatrixType(localMatrixStack_, domainEntity, rangeEntity);
      }

      LocalColumnObjectType localColumn( const DomainEntityType &colEntity ) const
      {
        return LocalColumnObjectType ( *this, colEntity );
      }

    public:
      void unitRow( const PetscInt localRow, const PetscScalar diag = 1.0 )
      {
        std::array< PetscInt, domainLocalBlockSize > rows;
        const PetscInt row = rangeMappers_.parallelIndex( localRow );
        for( unsigned int i=0, r = row * domainLocalBlockSize; i<domainLocalBlockSize; ++i, ++r )
        {
          rows[ i ] = r;
        }

        if( finalizedAlready() )
        {
          // set given row to a zero row with diagonal entry equal to diag
          ::Dune::Petsc::MatZeroRows( petscMatrix_, domainLocalBlockSize, rows.data(), diag );
        }
        else
        {
          // this only works for diag = 1
          assert( std::abs( diag - 1. ) < 1e-12 );
          unitRows_.reserve( unitRows_.size() + domainLocalBlockSize );
          for( const auto& r : rows )
          {
            unitRows_.push_back( r );
          }
        }
      }

      bool blockedMode() const { return blockedMode_; }

    protected:
      template< class PetscOp >
      void applyToBlock ( const PetscInt localRow, const PetscInt localCol, const MatrixBlockType& block, PetscOp op )
      {
#ifndef NDEBUG
        const PetscInt localCols = domainMappers_.ghostMapper().interiorSize() * domainLocalBlockSize;
        const PetscInt localRows = rangeMappers_.ghostMapper().interiorSize() * rangeLocalBlockSize;
        assert( localRow < localRows );
        assert( localCol < localCols );
#endif

        if( blockedMode_ )
        {
          // convert process local indices to global indices
          const PetscInt row = rangeMappers_.parallelIndex( localRow );
          const PetscInt col = rangeMappers_.parallelIndex( localCol );
          ::Dune::Petsc::MatSetValuesBlocked( petscMatrix_, 1, &row, 1, &col, block.data(), op );
        }
        else
        {
          // convert process local indices to global indices
          const PetscInt row = rangeMappers_.parallelIndex( localRow );
          const PetscInt col = rangeMappers_.parallelIndex( localCol );
          std::array< PetscInt, domainLocalBlockSize > rows;
          std::array< PetscInt, domainLocalBlockSize > cols;
          for( unsigned int i=0, r = row * domainLocalBlockSize, c = col * domainLocalBlockSize; i<domainLocalBlockSize; ++i, ++r, ++c )
          {
            rows[ i ] = r;
            cols[ i ] = c;
          }

          // set given row to a zero row with diagonal entry equal to diag
          ::Dune::Petsc::MatSetValues( petscMatrix_, domainLocalBlockSize, rows.data(), domainLocalBlockSize, cols.data(), block.data(), op );
        }
        setStatus( statAssembled );
      }

      template< class LocalBlock, class PetscOp >
      void applyToBlock ( const size_t row, const size_t col, const LocalBlock& block, PetscOp op )
      {
        assert( block.rows() == rangeLocalBlockSize );
        assert( block.cols() == domainLocalBlockSize );

        // copy to MatrixBlockType data structure suited to be inserted into Mat
        MatrixBlockType matBlock( block );
        applyToBlock( row, col, matBlock, op );
      }

      template< class LocalBlock >
      void setBlock ( const size_t row, const size_t col, const LocalBlock& block )
      {
#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statInsert );
#endif
        assert( row < std::numeric_limits< int > :: max() );
        assert( col < std::numeric_limits< int > :: max() );

        setStatus( statInsert );
        applyToBlock( static_cast< PetscInt > (row), static_cast< PetscInt > (col), block, INSERT_VALUES );
      }

      template< class LocalBlock >
      void addBlock ( const size_t row, const size_t col, const LocalBlock& block )
      {
#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statInsert );
#endif
        assert( row < std::numeric_limits< int > :: max() );
        assert( col < std::numeric_limits< int > :: max() );

        setStatus( statAdd );
        applyToBlock( static_cast< PetscInt > (row), static_cast< PetscInt > (col), block, ADD_VALUES );
      }

    protected:
      typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;

      // specialization for temporary local matrix, then copy of values is not needed
      template <class Operation>
      const PetscScalar* getValues( const unsigned int rSize, const unsigned int cSize,
                                    const TemporaryLocalMatrixType &localMat, const Operation&,
                                    const std::integral_constant< bool, false> nonscaled )
      {
        return localMat.data();
      }

      // specialization for temporary local matrix, then copy of values is not needed
      template <class LM, class Operation>
      const PetscScalar* getValues( const unsigned int rSize, const unsigned int cSize,
                                    const Assembly::Impl::LocalMatrixGetter< LM >& localMat, const Operation&,
                                    const std::integral_constant< bool, false> nonscaled )
      {
        return localMat.localMatrix().data();
      }

      // retrieve values for arbitrary local matrix
      template <class LocalMatrix, class Operation, bool T>
      const PetscScalar* getValues( const unsigned int rSize, const unsigned int cSize,
                                    const LocalMatrix &localMat, const Operation& operation,
                                    const std::integral_constant< bool, T> scaled )
      {
        std::vector< PetscScalar >& v = *(v_);
        v.resize( rSize * cSize );
        for( unsigned int i = 0, ic = 0 ; i< rSize; ++i )
        {
          for( unsigned int j =0; j< cSize; ++j, ++ic )
          {
            v[ ic ] = operation( localMat.get( i, j ) );
          }
        }
        return v.data();
      }

      template< class LocalMatrix, class Operation, class PetscOp, bool T >
      void applyLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity,
                              const LocalMatrix &localMat, const Operation& operation,
                              PetscOp petscOp,
                              const std::integral_constant<bool, T> scaled )
      {
        auto& rcTemp = *(rcTemp_);
        std::vector< PetscInt >& r = rcTemp.first;
        std::vector< PetscInt >& c = rcTemp.second;

        if( blockedMode_ )
        {
          setupIndicesBlocked( rangeMappers_,  rangeEntity,  r );
          setupIndicesBlocked( domainMappers_, domainEntity, c );

          // domainLocalBlockSize == rangeLocalBlockSize
          const unsigned int rSize = r.size() * domainLocalBlockSize ;
          const unsigned int cSize = c.size() * domainLocalBlockSize ;

          const PetscScalar* values = getValues( rSize, cSize, localMat, operation, scaled );
          ::Dune::Petsc::MatSetValuesBlocked( petscMatrix_, r.size(), r.data(), c.size(), c.data(), values, petscOp );
        }
        else
        {
          setupIndices( rangeMappers_,  rangeEntity,  r );
          setupIndices( domainMappers_, domainEntity, c );

          const unsigned int rSize = r.size();
          const unsigned int cSize = c.size();

          const PetscScalar* values = getValues( rSize, cSize, localMat, operation, scaled );
          ::Dune::Petsc::MatSetValues( petscMatrix_, rSize, r.data(), cSize, c.data(), values, petscOp );
        }
        setStatus( statAssembled );
      }

    public:
      template< class LocalMatrix >
      void addLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statAdd );
#endif
        setStatus( statAdd );

        auto operation = [] ( const PetscScalar& value ) -> PetscScalar { return value; };

        applyLocalMatrix( domainEntity, rangeEntity, localMat, operation, ADD_VALUES, std::integral_constant< bool, false>() );
      }

      template< class LocalMatrix, class Scalar >
      void addScaledLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat, const Scalar &s )
      {
#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statAdd );
#endif
        setStatus( statAdd );

        auto operation = [ &s ] ( const PetscScalar& value ) -> PetscScalar { return s * value; };

        applyLocalMatrix( domainEntity, rangeEntity, localMat, operation, ADD_VALUES, std::integral_constant< bool, true>() );
      }

      template< class LocalMatrix >
      void setLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const LocalMatrix &localMat )
      {
#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statInsert );
#endif
        setStatus( statInsert );

        auto operation = [] ( const PetscScalar& value ) -> PetscScalar { return value; };

        applyLocalMatrix( domainEntity, rangeEntity, localMat, operation, INSERT_VALUES, std::integral_constant< bool, false>() );
      }

      template< class LocalMatrix >
      void getLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, LocalMatrix &localMat ) const
      {
        // make sure matrix is in correct state before using
        finalizeAssembly();

#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statGet );
#endif
        setStatus( statGet );

        auto& rcTemp = *(rcTemp_);
        std::vector< PetscInt >& r = rcTemp.first;
        std::vector< PetscInt >& c = rcTemp.second;
        std::vector< PetscScalar >& v = *(v_);

        setupIndices( rangeMappers_, rangeEntity, r );
        setupIndices( domainMappers_, domainEntity, c );

        v.resize( r.size() * c.size() );
        ::Dune::Petsc::MatGetValues( petscMatrix_, r.size(), r.data(), c.size(), c.data(), v.data() );

        for( std::size_t i =0 ; i< r.size(); ++i )
          for( std::size_t j =0; j< c.size(); ++j )
            localMat.set( i, j, v[ i * c.size() +j ] );

        setStatus( statAssembled );
      }

      void scaleLocalMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity, const RangeFieldType &s )
      {
        // make sure matrix is in correct state before using
        finalizeAssembly();

#ifndef USE_SMP_PARALLEL
        assert( status_==statAssembled || status_==statGet );
#endif
        setStatus( statGet );

        auto& rcTemp = *(rcTemp_);
        std::vector< PetscInt >& r = rcTemp.first;
        std::vector< PetscInt >& c = rcTemp.second;
        std::vector< PetscScalar >& v = *(v_);

        setupIndices( rangeMappers_, rangeEntity, r );
        setupIndices( domainMappers_, domainEntity, c );

        const unsigned int rSize = r.size();
        const unsigned int cSize = c.size();
        const std::size_t size = rSize * cSize;

        v.resize( size );
        // get values
        ::Dune::Petsc::MatGetValues( petscMatrix_, r.size(), r.data(), c.size(), c.data(), v.data() );

        // scale values
        for( std::size_t i=0; i<size; ++i )
          v[ i ] *= s;

        // set values again
        ::Dune::Petsc::MatSetValues( petscMatrix_, rSize, r.data(), cSize, c.data(), v.data(), INSERT_VALUES );
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

      // return reference to PETSc matrix object
      Mat& exportMatrix () const
      {
        // make sure matrix is in correct state
        finalizeAssembly();
        return petscMatrix_;
      }

    private:
      PetscLinearOperator ();

      //! destructor deleting PETSc Mat object
      void destroy ()
      {
        if( status_ != statNothing )
        {
          ::Dune::Petsc::MatDestroy( &petscMatrix_ );
          status_ = statNothing ;
        }
        sequence_ = -1;
        unitRows_.clear();
        auxRows_.clear();
      }

      void setStatus (const Status &newstatus) const
      {
        // in case OpenMP is used simply avoid status check
#ifndef USE_SMP_PARALLEL
        status_ = newstatus;
#endif
      }

      template< class DFS, class Entity >
      static void setupIndices ( const PetscMappers< DFS > &mappers, const Entity &entity, std::vector< PetscInt > &indices )
      {
        NonBlockMapper< const typename PetscMappers< DFS >::ParallelMapperType, DFS::localBlockSize > nonBlockMapper( mappers.parallelMapper() );
        nonBlockMapper.map( entity, indices );
      }

      template< class DFS, class Entity >
      static void setupIndicesBlocked ( const PetscMappers< DFS > &mappers, const Entity &entity, std::vector< PetscInt > &indices )
      {
        mappers.parallelMapper().map( entity, indices );
      }

      /*
       * data fields
       */
      DomainMappersType domainMappers_;
      RangeMappersType  rangeMappers_;

      int sequence_ = -1;
      mutable Mat petscMatrix_;

      mutable LocalMatrixStackType localMatrixStack_;
      mutable Status status_;

      const bool viennaCL_;
      const bool blockedMode_;

      mutable std::unique_ptr< PetscDomainFunctionType > petscArg_;
      mutable std::unique_ptr< PetscRangeFunctionType  > petscDest_;

      mutable ThreadSafeValue< std::vector< PetscScalar > > v_;
      mutable ThreadSafeValue< std::pair< std::vector< PetscInt >, std::vector< PetscInt > > > rcTemp_;

      mutable std::vector< PetscInt > unitRows_;
      mutable std::vector< PetscInt > auxRows_;
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
      [[deprecated("Use TemporaryLocal Matrix and {add,set,get}LocalMatrix")]]
      LocalMatrix ( const PetscLinearOperatorType &petscLinOp,
                    const DomainSpaceType &domainSpace,
                    const RangeSpaceType &rangeSpace )
      : BaseType( domainSpace, rangeSpace ),
        petscLinearOperator_( petscLinOp )
      {}

      void finalize()
      {
      }

      void init ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
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

        status_ = statAssembled;
        petscLinearOperator_.setStatus(status_);
      }

      inline void add ( const int localRow, const int localCol, const RangeFieldType &value )
      {
        assert( status_==statAssembled || status_==statAdd );
        status_ = statAdd;
        petscLinearOperator_.setStatus(status_);
        ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) , value, ADD_VALUES );
      }

      inline void set(const int localRow, const int localCol, const RangeFieldType &value )
      {
        assert( status_==statAssembled || status_==statInsert );
        status_ = statInsert;
        petscLinearOperator_.setStatus(status_);
        ::Dune::Petsc::MatSetValue( petscMatrix(), globalRowIndex( localRow ), globalColIndex( localCol ) , value, INSERT_VALUES );
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

      inline void scale ( const RangeFieldType &factor )
      {
        const_cast< PetscLinearOperatorType& > (petscLinearOperator_).scaleLocalMatrix( this->domainEntity(), this->rangeEntity(), factor );
      }

    private:
      LocalMatrix ();

      Mat& petscMatrix () { return petscLinearOperator_.petscMatrix_; }
      const Mat& petscMatrix () const { return petscLinearOperator_.petscMatrix_; }

    public:
      int rows()    const { return rowIndices_.size(); }
      int columns() const { return colIndices_.size(); }

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
      mutable Status status_;
    };


  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCLINEAROPERATOR_HH
