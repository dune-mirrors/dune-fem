#ifndef DUNE_CACHE_TEST_HH
#define DUNE_CACHE_TEST_HH

#include "../cacheprovider.hh"


namespace Dune {
  namespace Fem {

  class CacheProvider_Test
  {
  public:
    virtual void run();

    void prismTest();
    void hexaTest();
    void tetraTest();
    void triangleTest();
    void quadTest();
  };

  } // end namespace Fem
} // end namespace Dune

#endif
