#include <config.h>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{
  namespace Fem
  {
    std::unique_ptr<MPIManager> MPIManager::instance_ = nullptr;
  }
}
