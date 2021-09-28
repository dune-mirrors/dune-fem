#ifndef DUNE_FEM_MPIMANAGER_HH
#define DUNE_FEM_MPIMANAGER_HH

#include <memory>

#include <dune/common/parallel/mpicommunication.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/fem/quadrature/caching/registry.hh>

#if HAVE_PETSC
#include <dune/fem/misc/petsc/petsccommon.hh>
#endif

#include <dune/fem/storage/singleton.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

namespace Dune
{

  namespace Fem
  {

    struct MPIManager
    {
      typedef Dune::CollectiveCommunication< MPIHelper::MPICommunicator >
        CollectiveCommunication;
    private:
      static MPIManager &instance ()
      {
        return Singleton< MPIManager > :: instance();
      }

      static bool mpiFinalized ()
      {
        bool finalized = false ;
#if HAVE_MPI
        // check that MPI was not already finalized
        {
          int wasFinalized = -1;
          MPI_Finalized( &wasFinalized );
          finalized = bool( wasFinalized );
        }
#endif // #if HAVE_MPI
        return finalized ;
      }

    public:
      //! destructor calling finalize if this has not been done
      ~MPIManager()
      {
        _finalize();
      }

      void _finalize()
      {
        if( ! mpiFinalized() )
        {
          if( petscWasInitializedHere_ )
            ::Dune::Petsc::finalize();

          // if MPI_Init was called here and finalize has not been
          // called yet, then this is the place to call it
          if( wasInitializedHere_ )
          {
#if HAVE_MPI
            MPI_Finalize();
#endif
          }
        }
      }

      static void finalize()
      {
        instance()._finalize();
      }

      static void initialize ( int &argc, char **&argv )
      {
        MPIHelper *&helper = instance().helper_;
        std::unique_ptr< CollectiveCommunication > &comm = instance().comm_;

        // initialize treading environment
        ThreadManager::initialize();

        // the following initialization overrides the MPI_Init in dune-common
        // to avoid a call to MPI_Finalize before all singletons have been deleted
#if HAVE_MPI
        int wasInitialized = -1;
        MPI_Initialized( &wasInitialized );
        if(!wasInitialized)
        {
#ifndef USE_SMP_PARALLEL
          // standard MPI_Init
          // call normal MPI_Init here to prevent MPIHelper to interfering
          // with MPI_Finalize one program exit which would cause failure
          {
            int is_initialized = MPI_Init(&argc, &argv);
            if( is_initialized != MPI_SUCCESS )
              DUNE_THROW(InvalidStateException,"MPI_Init failed!");
          }
#else     // threaded init
          {
            int provided;
            // use MPI_Init_thread for hybrid parallel programs
            int is_initialized = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided );

            if( is_initialized != MPI_SUCCESS )
              DUNE_THROW(InvalidStateException,"MPI_Init_thread failed!");

#if not defined NDEBUG && defined DUNE_DEVEL_MODE
            // for OpenMPI provided seems to be MPI_THREAD_SINGLE
            // but the bybrid version still works. On BlueGene systems
            // the MPI_THREAD_FUNNELED is really needed
            if( provided != MPI_THREAD_FUNNELED )
            {
              if( provided == MPI_THREAD_SINGLE )
                dwarn << "MPI thread support = single (instead of funneled)!" << std::endl;
              else
                dwarn << "WARNING: MPI thread support = " << provided << " != MPI_THREAD_FUNNELED " << MPI_THREAD_FUNNELED << std::endl;
            }
#endif    // end NDEBUG
          }
#endif    // end USE_SMP_PARALLEL
          instance().wasInitializedHere_ = true;

        } // end if(!wasInitialized)
#endif  // end HAVE_MPI

        // if already initialized, do nothing further
        if( helper && comm )
          return ;

        // this will just initialize the static variables inside MPIHelper but
        // not call MPI_Init again
        helper = &MPIHelper::instance( argc, argv );
        comm.reset( new CollectiveCommunication( helper->getCommunicator() ) );

#if HAVE_PETSC
        // initialize PETSc if present
        // returns true if PETSc was initialized during this call
        instance().petscWasInitializedHere_ =
          ::Dune::Petsc::initialize( rank() == 0, argc, argv );
#endif

        // initialize static variables of QuadratureStorageRegistry
        QuadratureStorageRegistry::initialize();
      }

      static const CollectiveCommunication &comm ()
      {
        const std::unique_ptr< CollectiveCommunication > &comm = instance().comm_;
        if( !comm )
          DUNE_THROW( InvalidStateException, "MPIManager has not been initialized." );
        return *comm;
      }

      static int rank ()
      {
        return comm().rank();
      }

      static int size ()
      {
        return comm().size();
      }

    private:
      MPIHelper *helper_ = nullptr;
      std::unique_ptr< CollectiveCommunication > comm_;
      bool wasInitializedHere_ = false ;
      bool petscWasInitializedHere_ = false ;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MPIMANAGER_HH
