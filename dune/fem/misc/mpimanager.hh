#ifndef DUNE_FEM_MPIMANAGER_HH
#define DUNE_FEM_MPIMANAGER_HH

#include <dune/common/mpihelper.hh>

namespace Dune
{

  class MPIManager
  {
    typedef MPIManager ThisType;
    
    MPIHelper *helper_;
    
    MPIManager ()
    : helper_( 0 )
    {}

    // prohibit copying and assignment
    MPIManager ( const ThisType & );
    ThisType &operator= ( const ThisType & );
    
    static MPIManager &instance ()
    {
      static MPIManager singleton;
      return singleton;
    }
    
  public:
    static void initialize ( int &argc, char **&argv )
    {
      instance().helper_ = &MPIHelper :: instance( argc, argv );
    }

    static MPIHelper &helper ()
    {
      MPIHelper *const helper = instance().helper_;
      if( helper == 0 )
        DUNE_THROW( InvalidStateException, "MPIManager has not been initialized." );
      return *helper;
    }

    static int rank ()
    {
      return helper().rank();
    }

    static int size ()
    {
      return helper().size();
    }
  };
  
}

#endif
