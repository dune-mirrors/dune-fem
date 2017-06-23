#ifndef DUNE_FEM_PETSCDOFMAPPING_HH
#define DUNE_FEM_PETSCDOFMAPPING_HH

#include <iostream>
#include <map>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/misc/mpimanager.hh>

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>

namespace Dune
{
  namespace Fem
  {


    /* =================================================
     * class PetscDofMappings
     *
     * This class implements 2 mappings what we need for PETSc integration.
     *
     * 1) The local slave mapping which maps local dune dof indices to the local indices the PETSc Vec. This mapping also
     * knows about slave dofs.
     *
     * 2) The global mapping, that maps local dune dof indices to global PETSc indices.
     *
     * =================================================
     */
    template <class Space>
    class PetscDofMappings
    {
      typedef PetscDofMappings ThisType;

    public:
      typedef Space                                    DiscreteFunctionSpaceType;
      typedef std::vector< int >                       DofMappingType;
      typedef DofMappingType::size_type                IndexType;
      typedef typename DofMappingType::value_type      DofType;

      typedef PetscInt  GlobalDofType ;
      typedef Dune::DynamicVector< GlobalDofType >           GlobalDofMappingType ;

      PetscDofMappings ( const DiscreteFunctionSpaceType *space )
      : numOwnedDofBlocks_( 0 ),
        numSlaveBlocks_( 0 ),
        processStartIndex_( 0 ),
        localSlaveMapping_(),
        globalDofMapping_(),
        sequence_( -1 )
      {
        // update dof mapping
        update(*space);
      }

      bool update (const DiscreteFunctionSpaceType &space)
      {
        auto &slaveDofs = space.slaveDofs();
        const int sequence = space.sequence();
        if( sequence_ != sequence )
        {
          initialize( space, slaveDofs );
          sequence_ = sequence ;
          return true ;
        }
        return false ;
      }


      GlobalDofType numOwnedDofBlocks () const { return numOwnedDofBlocks_; }

      GlobalDofType numSlaveBlocks () const { return numSlaveBlocks_; }

      GlobalDofType processStartIndex () const { return processStartIndex_; }

      size_t size () const { return localSlaveMapping_.size(); }

      const DofType& localSlaveMapping ( IndexType index ) const
      {
        return localSlaveMapping_[ index ];
      }

      const GlobalDofType& globalMapping ( const IndexType index ) const
      {
        return globalDofMapping_[ index ];
      }

      // is the dof with global DUNE index 'i' a slave dof?
      bool isSlave ( IndexType i ) const
      {
        return static_cast< int >( localSlaveMapping( i ) ) >= numOwnedDofBlocks();
      }


      ///////////////////////////////////////////////
      //  interface methods for PersistentIndexSet
      ///////////////////////////////////////////////
      template < class Stream >
      void write( const Stream& ) const {}

      template < class Stream >
      void read( const Stream& ) {}

      // void resize () { update (); }
      // bool compress() { return update (); }

      PetscDofMappings () = delete;
      PetscDofMappings ( const ThisType& ) = delete;
      PetscDofMappings& operator= ( const ThisType& ) = delete;

    private:
      template< typename SlaveDofsProvider >
      void initializeMappings ( const DiscreteFunctionSpaceType &space, const SlaveDofsProvider& slaveDofs )
      {
        // How the local slave mapping is build:
        // Let s_1 < ... < s_n be the slave dof indices (as given by the slaveDofs object) and let
        // d_1 < ... < d_m be the dof indices of dofs that are owned by this proces. Petsc thinks of slave dofs as
        // dofs in a PETSc Vec 'behind the array'. So the local slave mapping is now simply a vector with the following
        // components:
        //    d_1, d_2, ..., d_n, s_1, s_2, ..., s_n

        #ifndef NDEBUG
          int ownedDofBlocks = 0;
        #endif

        localSlaveMapping_.resize( space.blockMapper().size(), -1 );
        globalDofMapping_.resize( space.blockMapper().size(), 0 );

        std::fill( localSlaveMapping_.begin(), localSlaveMapping_.end(), -1 );
        std::fill( globalDofMapping_.begin(), globalDofMapping_.end(), 0 );

        // global dof index
        GlobalDofType index = processStartIndex_ ;
        for( int slave = 0, i = 0; slave < slaveDofs.size(); ++slave )
        {
          const int nextSlave = slaveDofs[ slave ];
          for(; i < nextSlave; ++i )
          {
            localSlaveMapping_[ i ] = i - slave;
            globalDofMapping_[ i ] = index++;

            #ifndef NDEBUG
              ++ownedDofBlocks;
            #endif
          }

          // omit the last slave
          if( static_cast< size_t >( i ) == localSlaveMapping_.size() )
            break;

          assert( static_cast< size_t >( i ) < localSlaveMapping_.size() );
          // skip the slave dof
          localSlaveMapping_[ i ] = numOwnedDofBlocks_ + slave;
          globalDofMapping_[ i ] = 0;
          ++i;
        }

        #ifndef NDEBUG
          checkDofMappingConsistency();
        #endif
        assert( numOwnedDofBlocks_ == ownedDofBlocks );

        typedef typename DiscreteFunctionSpaceType :: template ToNewDimRange< 1 > :: Type DofSpaceType ;
        typedef typename DiscreteFunctionSpaceType :: GridPartType  GridPartType;
        GridPartType& gridPart = const_cast< GridPartType& > (space.gridPart());

        if( space.continuous() )
        {
          typedef DofSpaceType GlobalDofSpaceType ;
          GlobalDofSpaceType dofSpace( gridPart );

          // store global dofs as a discrete function to use the already
          // built communication patterns to sum up the global dofs
          // which in this case simply sets the global numbers of the dofs
          // from the other ranks (we have to use the space's range field type)
          VectorDiscreteFunction< GlobalDofSpaceType, GlobalDofMappingType  >
            dofMappingFunc( "globalDofs", dofSpace, globalDofMapping_ );

          // do communication
          dofMappingFunc.communicate();
        }
        else
        {
          // for discontinuous solution we only need one dof per element --> FV space
          typedef FiniteVolumeSpace< typename DofSpaceType :: FunctionSpaceType,
                                     typename DofSpaceType :: GridPartType >  GlobalDofSpaceType ;
          GlobalDofSpaceType dofSpace( gridPart );

          // store global dofs as a discrete function to use the already
          // built communication patterns to sum up the global dofs
          // which in this case simply sets the global numbers of the dofs
          // from the other ranks (we have to use the space's range field type)
          VectorDiscreteFunction< GlobalDofSpaceType, GlobalDofMappingType  >
            dofMappingFunc( "globalDofs", dofSpace, globalDofMapping_ );

          // do communication
          dofMappingFunc.communicate();
        }
      }

      template< typename SlaveDofsProvider >
      void initialize ( const DiscreteFunctionSpaceType &space, const SlaveDofsProvider& slaveDofs )
      {
        numSlaveBlocks_    = slaveDofs.size() - 1;
        numOwnedDofBlocks_ = space.blockMapper().size() - numSlaveBlocks_;

        // start with index 0 (use unsigned long as buffers)
        unsigned long processStartIndex = 0;
        unsigned long numOwnedDofBlocks = numOwnedDofBlocks_;

        // initialize count start index for each process
        MPI_Exscan( &numOwnedDofBlocks, &processStartIndex, 1, MPI_UNSIGNED_LONG, MPI_SUM, PETSC_COMM_WORLD );

        // store my start index
        processStartIndex_ = processStartIndex ;

        initializeMappings( space, slaveDofs );

        #ifndef NDEBUG
          checkSlaveConsistency( slaveDofs );
        #endif
      }

      //////////////////////////////////////////
      // Methods for consistency checking, only used when NDEBUG is not set
      //////////////////////////////////////////
      template< typename SlavesType >
      void checkSlaveConsistency ( const SlavesType& slaves ) const
      {
        for( size_t i = 0; i < size(); ++i )
        {
          assert( isSlave( i ) == slaves.isSlave( i ) );
        }
      }

      void checkDofMappingConsistency () const
      {
        // Check if the dof mapping is bijective and valid. This piece of code does not strive to be
        // efficient...
        std::map< DofType, bool > buf;
        for( std::vector< int >::const_iterator it = localSlaveMapping_.begin(); it != localSlaveMapping_.end(); ++it )
        {
          if( *it < 0 )
          {
            std::cerr << "localSlaveMapping_ has not been initialized completely on rank " << MPIManager::rank() << std::endl;
            assert( false );
          }

          if( buf.find( *it ) != buf.end() )
          {
            std::cerr << *it << " was found more than once in localSlaveMapping_ (on rank " << MPIManager::rank() << ")\n";
            assert( false );
          }
          else
          {
            buf[ *it ] = true;
          }
        }
      }

      ////////////////////////////////
      // data fields
      ////////////////////////////////
      GlobalDofType   numOwnedDofBlocks_; // number of blocks where this proc is master
      GlobalDofType   numSlaveBlocks_;    // number of slave blocks
      GlobalDofType   processStartIndex_; // Start index of this process' portion of the PETSc Vec.
      DofMappingType  localSlaveMapping_; // local mapping (used for Discrete Function)
      GlobalDofMappingType globalDofMapping_; // global mapping (needed for matrices)
      int sequence_ ;                     // sequence number of adaptive grid
    };

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCDOFMAPPING_HH
