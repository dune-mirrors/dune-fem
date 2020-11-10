#include <config.h>

#include <dune/fem/storage/singleton.hh>

namespace Dune
{
  namespace Fem
  {
    namespace detail
    {
      typename SingletonStorage::StorageType* SingletonStorage::storage_ = nullptr;
    }

  } // namespace Fem

} // namespace Dune
