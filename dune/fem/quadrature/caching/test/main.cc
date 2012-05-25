#include <config.h>

#include <dune/fem/misc/suite.hh>

#include "refelem_test.hh"
#include "cache_test.hh"
#include "twist_test.hh"
#include "point_test.hh"

using namespace Dune;

int main() {

  Suite suite("Tests for caching");
  suite.addTest(new Fem::ReferenceElement_Test());
  suite.addTest(new Fem::TwistProvider_Test());
  suite.addTest(new Fem::PointProvider_Test());
  suite.addTest(new Fem::CacheProvider_Test());

  suite.run();
  suite.report();

}
