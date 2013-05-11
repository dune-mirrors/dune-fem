#ifndef DUNE_FEM_MPIMANAGER_HH
#define DUNE_FEM_MPIMANAGER_HH

#include <dune/common/parallel/mpicollectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>

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
      MPIManager ()
      : helper_( 0 ),
        comm_( 0 )
      {}

      ~MPIManager ()
      {
        delete comm_;

#if HAVE_PETSC
        ::Dune::Petsc::finalize();
#endif
      }

      // prohibit copying and assignment
      MPIManager ( const MPIManager & );
      MPIManager &operator= ( const MPIManager & );
      
      static MPIManager &instance ()
      {
        static MPIManager instance;
        return instance;
      }

    public:
      static void initialize ( int &argc, char **&argv )
      {
        // the following initalization is only enabled for 
        // MPI-thread parallel programs 
#if HAVE_MPI && MPI_2 
#ifdef USE_SMP_PARALLEL 
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
#endif // end NDEBUG 

#endif // end USE_SMP_PARALLEL
#endif // end HAVE_MPI && MPI_2       

        // call MPIHelper before instance(), this way MPIHelper seems to be deleted 
        // before MPIManager, which is needed because of the PETSc::finalize call 
        MPIHelper *helper = &MPIHelper::instance( argc, argv );
        CollectiveCommunication *comm = new CollectiveCommunication( helper->getCommunicator() ); 

        if( (instance().helper_ != 0) || (instance().comm_ != 0) )
          DUNE_THROW( InvalidStateException, "MPIManager has already been initialized." );

        // if not already called, this will call MPI_Init 
        instance().helper_ = helper ;
        instance().comm_   = comm ;

#if HAVE_PETSC
        ::Dune::Petsc::initialize( rank()==0, argc, argv );
#endif
      }

      static MPIHelper &helper () DUNE_DEPRECATED
      {
        MPIHelper *const helper = instance().helper_;
        if( helper == 0 )
          DUNE_THROW( InvalidStateException, "MPIManager has not been initialized." );
        return *helper;
      }

      static const CollectiveCommunication &comm ()
      {
        const CollectiveCommunication *const comm = instance().comm_;
        if( comm == 0 )
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
      MPIHelper *helper_;
      CollectiveCommunication *comm_;
    };

  } // namespace Fem

#if DUNE_FEM_COMPATIBILITY 
  using Fem::MPIManager;
#endif // DUNE_FEM_COMPATIBILITY 

} // namespace Dune

#endif // #ifndef DUNE_FEM_MPIMANAGER_HH
