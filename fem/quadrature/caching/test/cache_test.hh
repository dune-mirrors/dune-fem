#ifndef DUNE_CACHE_TEST_HH
#define DUNE_CACHE_TEST_HH

#include "../cacheprovider.hh"
#include <dune/fem/misc/test.hh>


namespace Dune {

  class CacheProvider_Test : public Test 
  {
  public:
    virtual void run();

    void hexaTest();
    void tetraTest();
    void triangleTest();
    void quadTest();
  };

} // end namespace Dune

#endif
