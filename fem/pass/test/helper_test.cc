#include <dune/config.h>

#include "helper_test.hh"
#include <dune/common/typetraits.hh>
#include <dune/quadrature/fixedorder.hh>

namespace Dune {
  void PassHelper_Test::run() {
    filterTest();
    tupleConverterTest();
    functorTest();
  }

  void PassHelper_Test::filterTest() {
    typedef Tuple<int, double, char, bool> TupleType;
    typedef TupleType::FirstPair PairType;
    typedef Selector<0, 2>::Base MySelector;
    typedef Selector<1>::Base SingleSelector;
    typedef Filter<PairType, MySelector>::ResultType MyFilteredType;
    typedef Filter<PairType, SingleSelector>::ResultType SingleFilteredType;

    TupleType t(0, 1.0, '3', true);
    MyFilteredType f = Filter<PairType, MySelector>::apply(t);

    _test(Element<0>::get(f) == Element<0>::get(t));
    _test(Element<1>::get(f) == Element<2>::get(t));

    SingleFilteredType fs = Filter<PairType, SingleSelector>::apply(t);
    _test(Element<0>::get(fs) == Element<1>::get(t));
  }

  void PassHelper_Test::tupleConverterTest() {
    typedef Lagrange_Fixture<0> Fix0;
    typedef Lagrange_Fixture<1> Fix1;
    typedef Lagrange_Fixture<0> Fix2;

    typedef DFAdapt<Fix0::DiscreteFunctionSpaceType> DF0;
    typedef DFAdapt<Fix1::DiscreteFunctionSpaceType> DF1;
    typedef DFAdapt<Fix2::DiscreteFunctionSpaceType> DF2;

    typedef DF0::LocalFunctionType LF0;
    typedef DF1::LocalFunctionType LF1;
    typedef DF2::LocalFunctionType LF2;

    typedef LF0::RangeType R0;
    typedef LF1::RangeType R1;
    typedef LF2::RangeType R2;

    typedef Tuple<DF0*, DF1*, DF2*> DFTupleType;
    typedef DFTupleType::FirstPair DFPairType;

    typedef LocalFunctionCreator<DFPairType>::ResultType LFTupleType;
    typedef Tuple<LF0, LF1, LF2>::FirstPair LFPairType;
    
    typedef RangeVectorCreator<LFPairType>::ResultType RTupleType;
    typedef Tuple<R0, R1, R2>::FirstPair RPairType;

    typedef Fix0::GridType GridType;

    GridType grid(gridFile_.c_str());
    Fix0 fix0(grid);
    Fix1 fix1(grid);
    Fix2 fix2(grid);

    DF0 df0("null", fix0.space());
    DF1 df1("eins", fix1.space());
    DF2 df2("zwei", fix2.space());

    DFTupleType dft(&df0, &df1, &df2);
    
    LFTupleType lft = LocalFunctionCreator<DFPairType>::apply(dft);
    //LocalFunctionCreator<DFPairType> lfc(dft);
    //   LFTupleType lft = lfc.evaluate();

    const bool sameLf = SameType<LFPairType, LFTupleType>::value;
    _test(sameLf == true);

    RTupleType rt = RangeVectorCreator<LFPairType>::apply(lft);
    //RangeVectorCreator<LFPairType> rc(lft);
    //RTupleType rt = rc.evaluate();

    const bool sameR = SameType<RPairType, RTupleType>::value;
    _test(sameR == true);
  }

  void PassHelper_Test::functorTest() {
    typedef Lagrange_Fixture<0> Fix0;
    typedef Lagrange_Fixture<1> Fix1;
    typedef Lagrange_Fixture<0> Fix2;

    typedef DFAdapt<Fix0::DiscreteFunctionSpaceType> DF0;
    typedef DFAdapt<Fix1::DiscreteFunctionSpaceType> DF1;
    typedef DFAdapt<Fix2::DiscreteFunctionSpaceType> DF2;

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
    
    typedef RangeVectorCreator<LFPairType>::ResultType RTupleType;
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

    RTupleType rt = RangeVectorCreator<LFPairType>::apply(lft);

    ForEachValuePair<DFTupleType, LFTupleType> forEachDFandLf(dft, lft);
    ForEachValuePair<LFTupleType, RTupleType> forEachLFandR(lft, rt);

    LocalFunctionSetter<Entity> setter(*it);      
    forEachDFandLf.apply(setter);

    LocalFunctionEvaluateLocal<Entity, DomainType> evaluator(*it, x);
    forEachLFandR.apply(evaluator);   
  }

  void PassHelper_Test::callerTest() {
    _fail("not implemented");
  }

} // end namespace Dune
