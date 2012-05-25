#ifndef DUNE_FEM_COMMHELPER_HH
#define DUNE_FEM_COMMHELPER_HH

#include <dune/common/collectivecommunication.hh>
#include <dune/common/mpicollectivecommunication.hh>
#include <dune/common/mpihelper.hh>

namespace Dune
{

  template< class CollectiveComm >
  struct CollectiveCommunicationHelper;

  template< class C >
  struct CollectiveCommunicationHelper< CollectiveCommunication< C > >
  {
    typedef CollectiveCommunication< C > CollectiveCommunicationType;

    static CollectiveCommunicationType defaultCommunication ()
    {
      return CollectiveCommunicationType();
    }
  };

#if HAVE_MPI
  template<>
  struct CollectiveCommunicationHelper< CollectiveCommunication< MPI_Comm > >
  {
    typedef CollectiveCommunication< MPI_Comm > CollectiveCommunicationType;

    static CollectiveCommunicationType defaultCommunication ()
    {
      return CollectiveCommunicationType( MPIHelper :: getCommunicator() );
    }
  };
#endif

  typedef CollectiveCommunication< MPIHelper :: MPICommunicator >
    DefaultCollectiveCommunicationType;
  
}

#endif
