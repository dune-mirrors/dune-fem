#include <config.h>

#include "cache_test.hh"
#include "twist_test.hh"
#include "point_test.hh"

using namespace Dune;
using namespace Fem;

int main()
{

  Fem::TwistProvider_Test().run();
  Fem::PointProvider_Test().run();
  Fem::CacheProvider_Test().run();

  return 0;
}
