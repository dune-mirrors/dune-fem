// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCGHOSTARRAYBUILDER_HH
#define DUNE_FEM_PETSCGHOSTARRAYBUILDER_HH

#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/misc/mpimanager.hh>

#if HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/misc/petsc/petscdofmappings.hh>
#include <dune/fem/misc/petsc/petscslavedofprovider.hh>


namespace Dune
{

  namespace Fem
  {

    /*
     * forward declarations
     */
    template< typename, typename > class PetscDofCommunicator;



    /* =================================================
     * class PetscGhostArrayBuilder
     *
     * This class is responsable to build the array of indices that indicates the position
     * of ghost dofs in the PETSc vec. Once an instance is initialized properly, one can
     * access the array using array()
     * =================================================
     */
    template< typename SlaveProvider, typename PDofMapping >
    class PetscGhostArrayBuilder
    {
      typedef PetscGhostArrayBuilder< SlaveProvider, PDofMapping > ThisType;

      typedef std::vector< std::pair< int, int > > StorageType;
      struct StorageSort
      {
        template< typename IT > bool operator () ( IT a, IT b ) { return a.first < b.first; }
      };

    public:
      typedef SlaveProvider SlaveProviderType;
      typedef typename SlaveProviderType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef PDofMapping PetscDofMappingType;
      typedef PetscInt value_type;
      static const int blockSize = DiscreteFunctionSpaceType::localBlockSize;


      explicit PetscGhostArrayBuilder ( SlaveProviderType& slaveDofs,
                                        PetscDofMappingType& petscDofMapping )
      : array_( 0 ),
        arrayNeedsUpdate_( true )
      {
        if( slaveDofs.space().gridPart().comm().size() > 1 ) // YaspGrid and lagrange polorder = 2 fails unless we skip this in serial code
        {
          try
          {
            PetscDofCommunicator< DiscreteFunctionSpaceType, ThisType > handle( slaveDofs.space(), petscDofMapping, *this );
            slaveDofs.space().gridPart().communicate( handle, DiscreteFunctionSpaceType::GridPartType::indexSetInterfaceType, ForwardCommunication );
          }
          catch( const Dune::Exception &e )
          {
            std::cerr << e << std::endl;
            std::cerr << "Exception thrown in: " << __FILE__ << " line: " << __LINE__ << std::endl;
            abort();
          }
        }

        // TODO: This is a corner case, but it might be possible that
        // everything works fine but this assertion fails
        //
        // the test is being removed because it does not work correctly with
        // dg spaces - needs to be modified first
        // assert( !hasDuplicatedGhosts());
      }

      ~PetscGhostArrayBuilder ()
      {
        delete[] array_;
      }

      void clear ()
      {
        vector_.clear();
        arrayNeedsUpdate_ = true;
      }

      bool hasDuplicatedGhosts () const
      {
        std::map< int, bool > map;
        for( StorageType::const_iterator it = vector_.begin(); it != vector_.end(); ++it )
        {
          if( map.find( it->second ) != map.end() )
            std::cerr << it->second << " found more than once on rk " << MPIManager::rank() << "\n";

          map[ it->second ] = true;
        }
        return ( map.size() < vector_.size() );
      }

      value_type* array ()
      {
        updateArray();
        return array_;
      }

      value_type* array () const
      {
        updateArray();
        return array_;
      }

      const value_type& operator[] ( size_t i ) const { return vector_[ i ].second; }

      value_type& operator[] ( size_t i )
      {
        arrayNeedsUpdate_ = true;
        return vector_[ i ].second;
      }

      void push_back ( int index, int value )
      {
        vector_.push_back( std::make_pair( index, value ) );
        arrayNeedsUpdate_ = true;
      }

      size_t size () const { return vector_.size(); }

    private:
      PetscGhostArrayBuilder ( const ThisType& );
      ThisType& operator= ( const ThisType & );
      PetscGhostArrayBuilder ();

      void updateArray () const
      {
        if( arrayNeedsUpdate_ )
        {
          delete[] array_;
          array_ = new value_type[ blockSize * size() ];
          std::sort( vector_.begin(), vector_.end(), StorageSort() );
          for( size_t block = 0; block < size(); ++block )
          {
            for( size_t inBlock = 0; inBlock < static_cast< size_t >( blockSize ); ++inBlock )
            {
              array_[ block*blockSize + inBlock ] = ( vector_[ block ].second ) * blockSize + inBlock;
            }
          }
          arrayNeedsUpdate_ = false;
        }
      }


      mutable StorageType vector_;
      mutable value_type *array_;
      mutable bool arrayNeedsUpdate_;
    };

    /* =================================================
     * class PetscDofCommunicator
     *
     * This is a helper class for PetscGhostArrayBuilder which is responsable for the
     * communication necessary to build the ghost array.
     * =================================================
     */
    template< typename DFSpace, typename PGhostArrayBuilder >
    class PetscDofCommunicator
    : public CommDataHandleIF< PetscDofCommunicator< DFSpace, PGhostArrayBuilder >, int >
    {
      typedef PetscDofCommunicator< DFSpace, PGhostArrayBuilder > ThisType;
      typedef CommDataHandleIF< PetscDofCommunicator< DFSpace, PGhostArrayBuilder >, int > BaseType;

    public:

      typedef PGhostArrayBuilder PetscGhostArrayBuilderType;
      typedef DFSpace DiscreteFunctionSpaceType;
      typedef typename BaseType::DataType DataType;
      typedef PetscSlaveDofProvider< DiscreteFunctionSpaceType > PetscSlaveDofProviderType ;
      typedef typename PetscSlaveDofProviderType :: PetscDofMappingType PetscDofMappingType;
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;

      /*
       * methods
       */
      PetscDofCommunicator ( const DiscreteFunctionSpaceType &space,
                             PetscDofMappingType &petscDofMapping,
                             PetscGhostArrayBuilderType& ghostArrayBuilder )
      : space_( space ),
        mapper_( space_.blockMapper() ),
        petscDofMapping_( petscDofMapping ),
        ghostArrayBuilder_( ghostArrayBuilder ),
        duplicates_()
      {}

      bool contains ( int dim, int codim ) const
      {
        return mapper_.contains( codim );
      }

      bool fixedsize ( int dim, int codim ) const
      {
        return false;
      }

      template< typename MessageBuffer, typename Entity >
      void gather ( MessageBuffer &buffer, const Entity &entity ) const
      {
        const int numDofs = mapper_.numEntityDofs( entity );
        int counter = 0;
        std::vector< int > slaveDofsWithIndex( 2*numDofs, -1 );

        std::vector< size_t > globalDofs;
        globalDofs.resize( numDofs );
        AssignFunctor< std::vector< size_t> >  functor( globalDofs );
        mapper_.mapEachEntityDof( entity, functor );

        for( int i = 0; i < numDofs; ++i )
        {
          size_t &dof = globalDofs[i];
          assert( size_t( dof ) < petscDofMapping_.size() );
          const int petscDof = petscDofMapping_.localSlaveMapping( dof );

          // does this process own the dof?
          if( petscDof < petscDofMapping_.numOwnedDofBlocks() )
          {
            const int petscGlobalDof = petscDofMapping_.processStartIndex() + petscDof;

            // write the entity dof index and the PETSc dof
            slaveDofsWithIndex[ counter++ ] = i;
            slaveDofsWithIndex[ counter++ ] = petscGlobalDof;
          }
        }
        assert( counter % 2 == 0 );

        // write the number of communicated slave dofs first
        buffer.write( counter );

        assert( size_t( counter ) <= slaveDofsWithIndex.size() );
        for( int i = 0; i < counter; ++i )
        {
          assert( slaveDofsWithIndex[ i ] != -1 );
          buffer.write( slaveDofsWithIndex[ i ] );
        }
      }

      template< typename MessageBuffer, typename EntityType >
      void scatter ( MessageBuffer &buffer, const EntityType &entity, size_t n )
      {
        bool isDuplicate = false;
        int idx = space_.indexSet().index(entity);
        if (duplicates_.find(idx) != duplicates_.end())
          isDuplicate = true;
        else
          duplicates_.insert(idx);

        int twiceNumSlaveDofs;
        buffer.read( twiceNumSlaveDofs );

        assert( twiceNumSlaveDofs % 2 == 0 );
        assert( static_cast< size_t >( twiceNumSlaveDofs ) + 1 == n );
        assert( twiceNumSlaveDofs + 1 == int( n ) );


        const int numDofs = mapper_.numEntityDofs( entity );
        std::vector< size_t > globalDofs;
        globalDofs.resize( numDofs );
        AssignFunctor< std::vector< size_t> >  functor( globalDofs );
        mapper_.mapEachEntityDof( entity, functor );

        for( size_t i = 0; i < size_t( twiceNumSlaveDofs/2 ); ++i )
        {
          int petscDof, index;
          buffer.read( index );

          size_t &globalDof = globalDofs[ index ];
          assert( petscDofMapping_.isSlave( globalDof ) );

          buffer.read( petscDof );
          if ( !isDuplicate )
            ghostArrayBuilder_.push_back( globalDof, petscDof );
        }
      }

      //! return local dof size to be communicated
      template< typename Entity >
      size_t size ( const Entity &entity ) const
      {
        size_t size = 0;

        const int numDofs = mapper_.numEntityDofs( entity );

        std::vector< size_t > globalDofs;
        globalDofs.resize( numDofs );
        AssignFunctor< std::vector< size_t> >  functor( globalDofs );
        mapper_.mapEachEntityDof( entity, functor );

        for( int i = 0; i < numDofs; ++i )
        {
          if( !petscDofMapping_.isSlave( globalDofs[i] ) )
            ++size;
        }
        assert( size == 0 || size == size_t( numDofs ) );

        return 2*size + 1;
      }

    private:
      PetscDofCommunicator ();
      PetscDofCommunicator ( const ThisType& );
      ThisType operator= ( const ThisType& );


      const DiscreteFunctionSpaceType &space_;
      const BlockMapperType &mapper_;
      PetscDofMappingType &petscDofMapping_;
      PetscGhostArrayBuilderType& ghostArrayBuilder_;
      std::set<int> duplicates_;

    };


  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCGHOSTARRAYBUILDER_HH
