#ifndef DUNE_FEM_MPIMANAGER_HH
#define DUNE_FEM_MPIMANAGER_HH

#include <dune/common/mpicollectivecommunication.hh>
#include <dune/common/mpihelper.hh>

namespace Dune
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
      MPIHelper *&helper = instance().helper_;
      CollectiveCommunication *&comm = instance().comm_;
#if HAVE_MPI && MPI_2 
#ifdef USE_SMP_PARALLEL 
      int provided;
      // use MPI_Init_thread for hybrid parallel programs 
      int is_initialized = MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided );

      if( is_initialized != MPI_SUCCESS ) 
        DUNE_THROW(InvalidStateException,"MPI_Init_thread failed!");

#if not defined NDEBUG && defined DUNE_DEVEL_MODE
      if( provided != MPI_THREAD_SINGLE ) 
        DUNE_THROW(InvalidStateException,"Wrong MPI mode initialized " << provided );
      /*
      if( provided != MPI_THREAD_FUNNELED )
      {
        if( provided == MPI_THREAD_SINGLE )
          dwarn << "MPI thread support = single (instead of funneled)!" << std::endl;
        else
          dwarn << "WARNING: MPI thread support = " << provided << " != MPI_THREAD_FUNNELED " << MPI_THREAD_FUNNELED << std::endl;
      }
      */

#endif // end NDEBUG 
#endif // end USE_SMP_PARALLEL
#endif // end HAVE_MPI && MPI_2       

      if( (helper != 0) || (comm != 0) )
        DUNE_THROW( InvalidStateException, "MPIManager has already been initialized." );

      helper = &MPIHelper::instance( argc, argv );
      comm = new CollectiveCommunication( helper->getCommunicator() );
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

}

#endif // #ifndef DUNE_FEM_MPIMANAGER_HH
