#include <dune/config.h>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/common/typetraits.hh>
#include "helper_test.hh"
//#include <dune/fem/quadrature/elementquadrature.hh>

namespace Dune {
  void PassHelper_Test::run() {
    selectorTest();
    filterTest();
    tupleConverterTest();
  }

  void PassHelper_Test::selectorTest() {
    typedef Selector<0, 2, 1>::Base Selector0;
    typedef Selector<-2, -3>::Base NegativeSelector;

    const int s0 = MaxIndex<Selector0>::value;
    _test(s0 == 2);

    const int ns = MaxIndex<NegativeSelector>::value;
    _test(ns == -1);

    const int tmp = MaxIndex<double>::value;
    _test(tmp == -1);
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
    typedef Lagrange_Fixture<GridType, 1> Fix0;
    typedef Lagrange_Fixture<GridType, 1> Fix1;
    typedef Lagrange_Fixture<GridType, 2> Fix2;

    typedef AdaptiveDiscreteFunction<Fix0::DiscreteFunctionSpaceType> DF0;
    typedef AdaptiveDiscreteFunction<Fix1::DiscreteFunctionSpaceType> DF1;
    typedef AdaptiveDiscreteFunction<Fix2::DiscreteFunctionSpaceType> DF2;

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
    
    typedef Creator<RangeTypeEvaluator, LFPairType>::ResultType RTupleType;
    typedef Tuple<R0, R1, R2>::FirstPair RPairType;

    typedef Fix0::GridType GridType;

    GridPtr<GridType> gridPtr(gridFile_);
    GridType& grid = *gridPtr;
    Fix0 fix0(grid);
    Fix1 fix1(grid);
    Fix2 fix2(grid);

    DF0 df0("null", fix0.space());
    DF1 df1("eins", fix1.space());
    DF2 df2("zwei", fix2.space());

    DFTupleType dft(&df0, &df1, &df2);
    
    LFTupleType lft(LocalFunctionCreator<DFPairType>::apply(dft));
    //LocalFunctionCreator<DFPairType> lfc(dft);
    //   LFTupleType lft = lfc.evaluate();

    const bool sameLf = is_same<LFPairType, LFTupleType>::value;
    _test(sameLf == true);

    RTupleType rt = Creator<RangeTypeEvaluator, LFPairType>::apply(lft);
    //RangeVectorCreator<LFPairType> rc(lft);
    //RTupleType rt = rc.evaluate();

    const bool sameR = is_same<RPairType, RTupleType>::value;
    _test(sameR == true);
  }


} // end namespace Dune
