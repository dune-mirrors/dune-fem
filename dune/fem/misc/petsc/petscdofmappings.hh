#ifndef DUNE_FEM_MISC_PETSC_PETSCDOFMAPPINGS_HH
#define DUNE_FEM_MISC_PETSC_PETSCDOFMAPPINGS_HH

#include <cstddef>

#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <vector>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/misc/mpimanager.hh>

#if HAVE_PETSC
#include <dune/fem/misc/petsc/petsccommon.hh>
#endif // #if HAVE_PETSC

namespace Dune
{

  namespace Fem
  {

    // PetscDofMappings
    // ----------------

    /*
     * implements 2 mappings needed for PETSc integration:
     * 1) The local slave mapping which maps local dune dof indices to the local
     *    indices the PETSc Vec. This mapping also knows about slave dofs.
     * 2) The global mapping, that maps local dune dof indices to global PETSc indices.
     */
#if HAVE_PETSC
    template< class Space, class GlobalDof = PetscInt >
    class PetscDofMappings
    {
      typedef PetscDofMappings ThisType;

    public:
      typedef Space DiscreteFunctionSpaceType;

    private:
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

      struct DataHandle;

    public:
      typedef int DofType;
      typedef std::vector< DofType > DofMappingType;
      typedef DofMappingType::size_type IndexType;

      typedef GlobalDof GlobalDofType;
      typedef std::vector< GlobalDofType > GlobalDofMappingType;

      PetscDofMappings ( const DiscreteFunctionSpaceType *space )
      {
        // update dof mapping
        update(*space);
      }

      PetscDofMappings ( const ThisType & ) = delete;
      PetscDofMappings ( ThisType && ) = delete;

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      bool update ( const DiscreteFunctionSpaceType &space )
      {
        if( sequence_ == space.sequence() )
          return false;

        initialize( space );
        sequence_ = space.sequence();
        return true;
      }

      GlobalDofType numOwnedDofBlocks () const { return numOwnedDofBlocks_; }
      GlobalDofType numSlaveBlocks () const { return numSlaveBlocks_; }

      GlobalDofType processStartIndex () const { return processStartIndex_; }

      IndexType size () const { return localSlaveMapping_.size(); }

      const DofType &localSlaveMapping ( IndexType index ) const { return localSlaveMapping_[ index ]; }
      const GlobalDofType &globalMapping ( const IndexType index ) const { return globalDofMapping_[ index ]; }

      // is the dof with global DUNE index 'i' a slave dof?
      bool isSlave ( IndexType i ) const { return (static_cast< int >( localSlaveMapping( i ) ) >= numOwnedDofBlocks()); }

      ///////////////////////////////////////////////
      //  interface methods for PersistentIndexSet
      ///////////////////////////////////////////////
      template < class Stream >
      void write( const Stream & ) const
      {}

      template < class Stream >
      void read( const Stream & )
      {}

      // void resize () { update (); }
      // bool compress() { return update (); }

    private:
      void initializeMappings ( const DiscreteFunctionSpaceType &space );

      void initialize ( const DiscreteFunctionSpaceType &space )
      {
        const auto &slaveDofs = space.slaveDofs();

        numSlaveBlocks_ = slaveDofs.size() - 1;
        numOwnedDofBlocks_ = space.blockMapper().size() - numSlaveBlocks_;

        // start with index 0 (use unsigned long as buffers)
        unsigned long processStartIndex = 0;
        unsigned long numOwnedDofBlocks = numOwnedDofBlocks_;

        // initialize count start index for each process
        MPI_Exscan( &numOwnedDofBlocks, &processStartIndex, 1, MPI_UNSIGNED_LONG, MPI_SUM, PETSC_COMM_WORLD );

        // store my start index
        processStartIndex_ = processStartIndex;

        initializeMappings( space );

        assert( checkSlaveConsistency( slaveDofs ) );
      }

      template< class SlaveDofs >
      bool checkSlaveConsistency ( const SlaveDofs &slaveDofs ) const
      {
        for( IndexType i = 0; i < size(); ++i )
        {
          if( isSlave( i ) != slaveDofs.isSlave( i ) )
            return false;
        }
        return true;
      }

      bool checkDofMappingConsistency () const
      {
        // Check if the dof mapping is bijective and valid. This piece of code does not strive to be
        // efficient...
        std::set< DofType > locals;
        for( const int local : localSlaveMapping_ )
        {
          if( local < 0 )
          {
            std::cerr << "localSlaveMapping_ has not been initialized completely on rank " << MPIManager::rank() << std::endl;
            return false;
          }

          if( !locals.insert( local ).second )
          {
            std::cerr << local << " was found more than once in localSlaveMapping_ (on rank " << MPIManager::rank() << ")" << std::endl;
            return false;
          }
        }
        return true;
      }

      GlobalDofType numOwnedDofBlocks_ = 0;   // number of blocks where this proc is master
      GlobalDofType numSlaveBlocks_ = 0;      // number of slave blocks
      GlobalDofType processStartIndex_ = 0;   // Start index of this process' portion of the PETSc Vec.
      DofMappingType localSlaveMapping_;      // local mapping (used for Discrete Function)
      GlobalDofMappingType globalDofMapping_; // global mapping (needed for matrices)
      int sequence_ = -1;                     // sequence number of adaptive grid
    };



    // PetscDofMappings::DataHandle
    // ----------------------------

    template< class Space, class GlobalDof >
    struct PetscDofMappings< Space, GlobalDof >::DataHandle
      : public CommDataHandleIF< DataHandle, GlobalDof >
    {
      DataHandle ( const BlockMapperType &blockMapper, GlobalDofMappingType &globalDofMapping )
        : blockMapper_( blockMapper ), globalDofMapping_( globalDofMapping )
      {}

      bool contains ( int dim, int codim ) const { return blockMapper_.contains( codim ); }
      bool fixedsize ( int dim, int codim ) const { return blockMapper_.fixedDataSize( codim ); }

      template< class Buffer, class Entity >
      void gather ( Buffer &buffer, const Entity &entity ) const
      {
        blockMapper_.mapEachEntityDof( entity, [ this, &buffer ] ( int, auto i ) { buffer.write( globalDofMapping_[ i ] ); } );
      }

      template< class Buffer, class Entity >
      void scatter ( Buffer &buffer, const Entity &entity, std::size_t n )
      {
        assert( n == blockMapper_.numEntityDofs( entity ) );
        blockMapper_.mapEachEntityDof( entity, [ this, &buffer ] ( int, auto i ) {
            GlobalDofType index;
            buffer.read( index );
            globalDofMapping_[ i ] = std::min( globalDofMapping_[ i ], index );
          } );
      }

      template< class Entity >
      std::size_t size ( const Entity &entity ) const
      {
        return blockMapper_.numEntityDofs( entity );
      }

    private:
      const BlockMapperType &blockMapper_;
      GlobalDofMappingType &globalDofMapping_;
    };



    template< class Space, class GlobalDof >
    inline void PetscDofMappings< Space, GlobalDof >::initializeMappings ( const DiscreteFunctionSpaceType &space )
    {
      // How the local slave mapping is build:
      // Let s_1 < ... < s_n be the slave dof indices (as given by the slaveDofs object) and let
      // d_1 < ... < d_m be the dof indices of dofs that are owned by this proces. Petsc thinks of slave dofs as
      // dofs in a PETSc Vec 'behind the array'. So the local slave mapping is now simply a vector with the following
      // components:
      //    d_1, d_2, ..., d_n, s_1, s_2, ..., s_n

      const auto &slaveDofs = space.slaveDofs();
      const BlockMapperType &blockMapper = space.blockMapper();

      localSlaveMapping_.resize( blockMapper.size() );
      globalDofMapping_.resize( blockMapper.size() );

      std::fill( localSlaveMapping_.begin(), localSlaveMapping_.end(), -1 );
      std::fill( globalDofMapping_.begin(), globalDofMapping_.end(), std::numeric_limits< GlobalDofType >::max() );

      // global dof index
      GlobalDofType globalIndex = processStartIndex_;
      IndexType index = 0;
      for( int slave = 0; slave < slaveDofs.size(); ++slave )
      {
        const IndexType nextSlave = slaveDofs[ slave ];
        for( ; index < nextSlave; ++index )
        {
          localSlaveMapping_[ index ] = index - slave;
          globalDofMapping_[ index ] = globalIndex++;
        }

        // omit the last slave
        if( index == localSlaveMapping_.size() )
          break;
        assert( index < localSlaveMapping_.size() );

        // assign the slave DoF a ghost index
        localSlaveMapping_[ index++ ] = numOwnedDofBlocks_ + slave;
      }

      assert( index == localSlaveMapping_.size() );
      assert( globalIndex == processStartIndex_ + numOwnedDofBlocks_ );
      assert( checkDofMappingConsistency() );

      DataHandle dataHandle( blockMapper, globalDofMapping_ );
      space.gridPart().communicate( dataHandle, GridPartType::indexSetInterfaceType, ForwardCommunication );
    }
#endif // #if HAVE_PETSC

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MISC_PETSC_PETSCDOFMAPPINGS_HH
