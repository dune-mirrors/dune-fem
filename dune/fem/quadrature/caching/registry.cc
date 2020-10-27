#include <config.h>
#include <dune/fem/quadrature/caching/registry.hh>

namespace Dune
{
  namespace Fem
  {
    QuadratureStorageRegistry::StorageListType QuadratureStorageRegistry::storageList_;
    QuadratureStorageRegistry::QuadratureInfoListType QuadratureStorageRegistry::quadratureInfoList_;
  }
}
