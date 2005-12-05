#include "../../misc/suite.hh"

#include "refelem_test.hh"

using namespace Dune;

int main() {

  Suite suite("Tests for caching");
  suite.addTest(new ReferenceElement_Test());

  suite.run();
  suite.report();

}
