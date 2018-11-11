#ifndef DUNE_FEM_MPIMANAGER_HH
#define DUNE_FEM_MPIMANAGER_HH

#include <memory>

#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/visibility.hh>

#include <dune/fem/quadrature/caching/registry.hh>

#if HAVE_PETSC
#include <dune/fem/misc/petsc/petsccommon.hh>
#endif

namespace Dune
{

  namespace Fem
  {

    struct MPIManager
    {
      typedef Dune::CollectiveCommunication< MPIHelper::MPICommunicator >
        CollectiveCommunication;

    private:
      DUNE_EXPORT static MPIManager &instance ()
      {
        static MPIManager instance;
        return instance;
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

#if HAVE_PETSC
      struct PETSc
      {
        ~PETSc()
        {
          if( ! mpiFinalized() )
          {
            ::Dune::Petsc::finalize();
          }
        }

        DUNE_EXPORT static void initialize( const bool verbose, int &argc, char **&argv )
        {
          // needed for later calling Petsc::finalize to the right time
          static PETSc petsc;
          ::Dune::Petsc::initialize( verbose, argc, argv );
        }
      };
#endif // #if HAVE_PETSC

    public:
      static void initialize ( int &argc, char **&argv )
      {
        MPIHelper *&helper = instance().helper_;
        std::unique_ptr< CollectiveCommunication > &comm = instance().comm_;

        // the following initalization is only enabled for
        // MPI-thread parallel programs
#if HAVE_MPI
#ifdef USE_SMP_PARALLEL
        int wasInitialized = -1;
        MPI_Initialized( &wasInitialized );
        if(!wasInitialized)
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
#endif  // end NDEBUG
        } // end if(!wasInitialized)
#endif  // end USE_SMP_PARALLEL
#endif  // end HAVE_MPI

        // if already initialized, do nothing further
        if( helper && comm )
          return ;

        // if not already called, this will call MPI_Init
        helper = &MPIHelper::instance( argc, argv );
        comm.reset( new CollectiveCommunication( helper->getCommunicator() ) );

#if HAVE_PETSC
        // initialize PETSc if pressent
        PETSc::initialize( rank() == 0, argc, argv );
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
        //int wasFinalized = -1;
        //MPI_Finalized( &wasFinalized );
        //if( wasFinalized)
        //{
        //  assert(false);
        //  std::abort();
        //}

        return comm().rank();
      }

      static int size ()
      {
        return comm().size();
      }

    private:
      MPIHelper *helper_ = nullptr;
      std::unique_ptr< CollectiveCommunication > comm_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_MPIMANAGER_HH
