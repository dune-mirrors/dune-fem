#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/quadrature/caching/registry.hh>

namespace Dune
{
  namespace Fem
  {
    void MPIManager::initialize ( int &argc, char **&argv )
    {
      MPIHelper *&helper = instance().helper_;
      std::unique_ptr< CollectiveCommunication > &comm = instance().comm_;

      // initialize treading environment
      // MPIManager::initialize();

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
  }
}
