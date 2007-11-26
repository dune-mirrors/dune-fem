#include <config.h>

#include <dune/fem/misc/suite.hh>

#include "refelem_test.hh"
#include "cache_test.hh"
#include "twist_test.hh"
#include "point_test.hh"

using namespace Dune;

int main() {

  Suite suite("Tests for caching");
  suite.addTest(new ReferenceElement_Test());
  suite.addTest(new TwistProvider_Test());
  suite.addTest(new PointProvider_Test());
  suite.addTest(new CacheProvider_Test());

  suite.run();
  suite.report();

}
