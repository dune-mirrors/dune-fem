#include <dune/config.h>

#include "pass_test.hh"
#include "passstub.hh"
#include "dgstub.hh"
#include "lagrange_fixture.hh"

namespace Dune {

  void Pass_Test::run() {
    dummyTest();
    dgTest();
  }

  void Pass_Test::dummyTest() {
    typedef PassStubTraits::DestinationType GlobalArgumentType;
    typedef PassStub<StartPass<GlobalArgumentType> > Pass1Type;
    typedef PassStub<Pass1Type> Pass2Type;

    StartPass<GlobalArgumentType> sp;
    Pass1Type p1(sp);
    Pass2Type p2(p1);

    _fail("Not implemented yet");
  }

  void Pass_Test::dgTest() {
    typedef DGStubTraits::DestinationType GlobalArgumentType;
    typedef DGStubTraits::DestinationType DestinationType;
    //typedef DGStubTraits::SpaceType SpaceType;
    typedef StartPass<GlobalArgumentType> Pass0Type;
    typedef LocalDGPass<ProblemStub, DGStubTraits, Pass0Type> Pass1Type;
    typedef LocalDGPass<ProblemStub, DGStubTraits, Pass1Type> Pass2Type;
    typedef Lagrange_Fixture<1> Fix;
    typedef Fix::GridType GridType;
    typedef Fix::DiscreteFunctionSpaceType SpaceType;

    GridType grid(gridFile_.c_str());
    Fix fix(grid);
    SpaceType& spc = fix.space();
    ProblemStub ps;

    Pass0Type p0;
    Pass1Type p1(ps, p0, spc);
    Pass2Type p2(ps, p1, spc);

    GlobalArgumentType arg("Arg", spc);
    DestinationType dest("Dest", spc);

    p2(arg, dest);

    _fail("Not implemented yet");

  }
} // end namespace Dune
