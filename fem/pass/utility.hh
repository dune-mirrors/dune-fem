#ifndef DUNE_PASS_UTILITY_HH
#define DUNE_PASS_UTILITY_HH

#include <dune/common/tuples.hh>

namespace Dune {

  /** @ addtogroup Common
   *
   * @{
   */

  /**
   * @file
   * @brief Contain utility classes which can be used with tuples.
   */

  template <class T1,class T2>
  struct CombinePairs {
  };
  
  /*********************************************/
template <class TupType>
struct TupleToPair {
  typedef typename TupleToPair<typename TupType::Type2>::ReturnType NextType;
  typedef Pair<typename TupType::Type1,NextType> ReturnType;
  static ReturnType convert(TupType tup) {
    NextType next = TupleToPair<typename TupType::Type2>::convert(tup.second());
    return ReturnType(tup.first(),next);
  }
};
template <class T>
struct TupleToPair<Pair<T,Nil> > {
  typedef Pair<T,Nil> ReturnType;
  static ReturnType convert(Pair<T,Nil> tup) {
    return ReturnType(tup.first(),nullType());
  }
};
/*
template <class T1,class T2> 
struct AddToTuple {
  typedef Pair<T1,Pair<T2,Nil> > ResultType;
  static ReturnType make(T1 t2,T2 t2) {
    return ResultType(makePair(t1,makePair(t2,Nil()));
  }
};
template <class T1> 
struct AddToTuple<T1,Nil> {
  typedef Pair<T1,Nil> ResultType;
  static ReturnType make(T1 t1) {
    return ResultType(makePair(t1,Nil());
  }
};
template <class T1,class T2,class T3,class T4,class T5,class T6,class T7,
	  class T8,class T9> 
struct AddToTuple<Tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> {
  typedef Tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> TupleType;
  typedef TupleToPair<TupleType>::ReturnType ResultType;
  static ReturnType make(TupleType t) {
    return TupleToPair<TupleType>::apply(t);
  }
};
template <class T11,class T12,class T21,class T22> 
struct AddToTuple<Pair<T11,T12>,Pair<T21,T22> > {
  typedef Pair<T11,T12> P1;
  typedef Pair<T21,T22> P2;
  typedef CombinePairs<P1,P2>::ReturnType ResultType;
  static ReturnType make(P1 p1,P2 p2) {
    return CombinePairs<P1,P2>::apply(p1,p2);
  }
};
*/
}

#endif
