#include <config.h>

#include "refelem_test.hh"
#include "cache_test.hh"
#include "twist_test.hh"
#include "point_test.hh"

using namespace Dune;
using namespace Fem;

int main() {

  //Suite suite("Tests for caching");
  Fem::ReferenceElement_Test().run();
  Fem::TwistProvider_Test().run();
  Fem::PointProvider_Test().run();
  Fem::CacheProvider_Test().run();

  return 0;
}
