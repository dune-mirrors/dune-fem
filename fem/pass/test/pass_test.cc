#include <dune/config.h>

#include "pass_test.hh"
#include "passstub.hh"
#include "dgstub.hh"
#include "lagrange_fixture.hh"

#include <dune/fem/discretefunction/adaptivefunction.hh>

namespace Dune {

  void Pass_Test::run() {
    functorTest();
    dummyTest();
    dgTest();
  }

  void Pass_Test::functorTest() {
    typedef Lagrange_Fixture<0> Fix0;
    typedef Lagrange_Fixture<1> Fix1;
    typedef Lagrange_Fixture<1> Fix2;

    typedef DFAdapt<Fix0::DiscreteFunctionSpaceType> DF0;
    typedef DFAdapt<Fix1::DiscreteFunctionSpaceType> DF1;
    typedef DFAdapt<Fix2::DiscreteFunctionSpaceType> DF2;

    //typedef AdaptiveDiscreteFunction<Fix0::DiscreteFunctionSpaceType> DF0;
    //typedef AdaptiveDiscreteFunction<Fix1::DiscreteFunctionSpaceType> DF1;
    //typedef AdaptiveDiscreteFunction<Fix2::DiscreteFunctionSpaceType> DF2;

    typedef DF0::LocalFunctionType LF0;
    typedef DF1::LocalFunctionType LF1;
    typedef DF2::LocalFunctionType LF2;

    typedef LF0::RangeType R0;
    typedef LF1::RangeType R1;
    typedef LF2::RangeType R2;

    typedef LF0::DomainType DomainType;

    typedef Tuple<DF0*, DF1*, DF2*> DFTupleType;
    typedef DFTupleType::FirstPair DFPairType;

    typedef LocalFunctionCreator<DFPairType>::ResultType LFTupleType;
    typedef Tuple<LF0, LF1, LF2>::FirstPair LFPairType;
    
    typedef Creator<RangeTypeEvaluator, LFPairType>::ResultType RTupleType;
    typedef Tuple<R0, R1, R2>::FirstPair RPairType;

    typedef Fix0::GridType GridType;
    typedef GridType::Codim<0>::LeafIterator LeafIterator;
    typedef LeafIterator::Entity Entity;

    GridType grid(gridFile_.c_str());
    Fix0 fix0(grid);
    Fix1 fix1(grid);
    Fix2 fix2(grid);

    DF0 df0("null", fix0.space());
    DF1 df1("eins", fix1.space());
    DF2 df2("zwei", fix2.space());

    LeafIterator it = grid.leafbegin<0>();
    DomainType x(0.2);

    DF0::LocalFunctionType lf = df0.localFunction(*it);
     
    DFTupleType dft(&df0, &df1, &df2);
   
    LFTupleType lft = LocalFunctionCreator<DFPairType>::apply(dft);

    //LFTupleType lft2 = LocalFunctionCreator2::apply(dft);

    RTupleType rt = Creator<RangeTypeEvaluator, LFPairType>::apply(lft);

    ForEachValuePair<DFTupleType, LFTupleType> forEachDFandLf(dft, lft);
    ForEachValuePair<LFTupleType, RTupleType> forEachLFandR(lft, rt);

    LocalFunctionSetter<Entity> setter(*it);      
    forEachDFandLf.apply(setter);

    LocalFunctionEvaluateLocal<Entity, DomainType> evaluator(*it, x);
    forEachLFandR.apply(evaluator);
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
    typedef LocalDGPass<ProblemStub, Pass0Type> Pass1Type;
    typedef LocalDGPass<ProblemStub, Pass1Type> Pass2Type;
    typedef Lagrange_Fixture<1> Fix;
    typedef Fix::GridType GridType;
    typedef Fix::DiscreteFunctionSpaceType SpaceType;
    typedef BoundaryManager<Fix::DiscreteFunctionSpaceType> BCManager;

    GridType grid(gridFile_.c_str());
    Fix fix(grid);
    SpaceType& spc = fix.space();
    ProblemStub ps;
    BCManager bcManager;

    Pass0Type p0;
    Pass1Type p1(ps, p0, spc, bcManager);
    Pass2Type p2(ps, p1, spc, bcManager);

    GlobalArgumentType arg("Arg", spc);
    DestinationType dest("Dest", spc);

    p2(arg, dest);

    _fail("Not implemented yet");

  }
} // end namespace Dune
